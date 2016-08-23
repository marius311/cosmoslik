from collections import OrderedDict
import copy, pkgutil, inspect, sys, os
from numpy import inf, nan, hstack, transpose
import imp, hashlib, time
import inspect
import argparse
from functools import reduce
from uuid import uuid4

__all__ = ['load_script','Slik','SlikFunction',
           'SlikDict','SlikPlugin','SlikSampler','param','param_shortcut',
           'get_plugin','get_all_plugins',
           'lsum','lsumk','all_kw','run_chain','SlikMain', 'arguments']           

""" Loaded datafiles will reside in this empty module. """
datafile_module = 'cosmoslik.scripts'
sys.modules[datafile_module]=imp.new_module(datafile_module)


def load_script(script):   
    """ 
    Read a CosmoSlik script from a file and return the :class:`SlikPlugin`
    object.
    """
    if isinstance(script,str):
        script_module = uuid4().hex
        modname='%s.%s'%(datafile_module,script_module)
        mod = imp.load_source(modname,script)
        sys.modules[modname]=mod
        
        plugins = [v for k,v in list(vars(mod).items()) if (isinstance(v,type) and 
                                                      issubclass(v,SlikPlugin) and 
                                                      v!=SlikPlugin)]
        mains = [v for v in plugins if hasattr(v,'_slik_main')]
        
        if len(mains)>=2: raise ValueError("Multiple SlikPlugins in '%s' are marked with @SlikMain. CosmoSlike doesn't know which one to run."%script)
        elif len(mains)==1: main=mains[0]
        elif len(plugins)>=2: raise ValueError("Multiple SlikPlugins were found in '%s' but none are marked with @SlikMain. CosmoSlik doesn't know which one to run."%script)
        elif len(plugins)==1: main=plugins[0]
        else: raise ValueError("No SlikPlugins were found in '%s'"%script)
    elif issubclass(script,SlikPlugin):
        main = script
    else:
        raise ValueError("`script` argument should be filename or SlikPlugin class")
    
    argspec = inspect.getargspec(main.__init__)
    class NoDefault: pass
    parser = argparse.ArgumentParser(prog="cosmoslik %s"%script)
    args = argspec.args[1:]
    defaults = [NoDefault]*(len(args) - len(argspec.defaults or [])) + list(argspec.defaults or [])
    for name, default in zip(args,defaults):
        if default == NoDefault:
            parser.add_argument(name)
        else:
            if default is True:
                parser.add_argument("--no-"+name, dest=name, action="store_false")
            elif default is False:
                parser.add_argument("--"+name, action="store_true")
            elif isinstance(default,list):
                parser.add_argument("--"+name, default=default, nargs=('+' if len(default)>0 else '*'), 
                                    type=(type(default[0]) if len(default)>0 else None),
                                    help="default: "+str(default))
            else:
                parser.add_argument("--"+name, default=default, type=type(default), help="default: "+str(default))
                
    
    return parser, main
    # return Slik(main(*args,**kwargs))
    
    
    
    


class Slik(object):
    
    def __init__(self,params):
        self.params = params
        self._sampled = params.find_sampled()
        self.sampler = self.params.sampler
        del self.params.sampler
        
        
    def get_sampled(self):
        """
        Get all the sampled parameters. 
        """
        return self._sampled
        
    def get_start(self):
        """
        Get the starting point as a dictionary which can immediately be passed
        method:`evaluate`::
        
            slik = Slik(...)
            slik.evaluate(**slik.get_start())
            
        """
        return {k:v.start for k,v in self.get_sampled().items()}
        
    def evaluate(self,*args,**kwargs):
        """
        evaluate(*args) - args is a vector 
        containing the values of the parameters in the order given by 
        
        or evluate(**kwagrs) - kwargs is a dictionary
        containing a mapping of all the parameters to set
        
        returns (likelihood, params)
        """
        
        
        params = self.params.deepcopy()
        
        if len(args)>0 and len(kwargs)>0: raise ValueError("Expected either *args or **kwargs but not both.")
        
        if len(args)==0:
            missing = [k for k in self.get_sampled() if k not in kwargs]
            if len(missing)>0:
                raise ValueError("Missing the following parameters: %s"%missing)
            for k,v in list(kwargs.items()): params[k]=v
        elif len(args)!=len(self.get_sampled()):
            raise ValueError("Expected %i parameters, only got %i."%(len(self.get_sampled()),len(args)))
        else:
            for k,v in zip(list(self.get_sampled().keys()),args): params[k]=v
            
        return params(), params
            
            
    def sample(self):
        return self.sampler.sample(self.evaluate)
    
    


def SlikMain(cls):
    """
    Class decorator which marks a plugin as being the main one when running a
    script from the command line.
    
    Example::
    
        # somefile.py
        
        @SlikMain  #run this plugin when I run `$ cosmoslik somefile.py`
        class plugin1(SlikPlugin):
            ...
        
        class plugin2(SlikPlugin):
            ...
            
    """
    if not issubclass(cls,SlikPlugin): raise ValueError("SlikMain used on a class which is not a SlikPlugin.")
    cls._slik_main = True
    return cls
        
        
def SlikFunction(func):
    func._slik_function = True
    return func


class SlikDict(dict):
    
    def __init__(self,*args,**kwargs):
        super().__init__(*args, **kwargs)
        self.__dict__ = self
        
    def __setitem__(self,k,v):
        if isinstance(k,str):
            setattr(reduce(getattr, k.split('.')[:-1], self), k.split('.')[-1], v)
        else:
            raise ValueError("Parameter key must be string, not %s"%type(k))

    def deepcopy(self):
        cself = copy.copy(self)
        for k,v in list(vars(self).items()):
            if isinstance(v,SlikDict): setattr(cself,k,v.deepcopy())
        cself.update(vars(cself))
        cself.__dict__ = cself
        return cself
    
    def __getitem__(self,k):
        if isinstance(k,str):
            return reduce(getattr,k.split('.'),self)
        else:
            raise ValueError("Parameter key must be string, not %s"%type(k))
    
    def get(self,k,default=None):
        if isinstance(k,str):
            return reduce(lambda obj,k: getattr(obj,k,default),k.split('.'),self)
        else:
            raise ValueError("Parameter key must be string, not %s"%type(k))

    def find_sampled(self):
        
        def _find_sampled(self,root):
            all_sampled = {}
            for k,v in self.__dict__.items(): 
                if isinstance(v,SlikDict): 
                    all_sampled.update(_find_sampled(v,root=root+[k]))
                elif isinstance(v,param):
                    all_sampled[('.'.join(root+[k]))]=v 
            return all_sampled
        
        return OrderedDict(sorted([(k,v) for k,v in _find_sampled(self,[]).items()]))


class param(object):
    def __init__(self,**kwargs):
        self.__dict__.update(kwargs)


class sample(object):
    def __init__(self,x,lnl,extra):
        self.x = x
        self.lnl = lnl
        self.extra = extra



def run_chain(main,nchains=1,pool=None,args=None,kwargs=None):
    """
    Run a CosmoSlik chain, or if `nchains`!=1 run a set of chains in parallel.

    Non-trivial parallelization, e.g. using MPI or with communication amongst
    chains, is handled by each cosmoslik.Sampler module. See individual
    documentation. 


    Args:
        main: the class object for your :class:`SlikPlugin` (i.e. ``myplugin``,
            not ``myplugin()``)
        pool: any worker pool which has a ``pool.map`` function. Defaults to
            ``multiprocessing.Pool(nchains)``
        nchains(int): the number of chains to run in parallel using ``pool.map``
        *args: args to pass to `main`
        **kwargs: kwargs to pass to `main`


    Returns:
        A :class:`chain.Chain` instance if nchains=1, 
        otherwise a :class:`chain.Chains` instance

    """
    from .chains import Chain, Chains
    from multiprocessing.pool import Pool

    if args is None: args=[]
    if kwargs is None: kwargs={}
    if isinstance(main,tuple): 
        main,args,kwargs = main

    if nchains==1:
        slik = Slik(main(*args,**kwargs))
        return Chain(dict(list(zip(hstack(['weight','lnl',list(slik.get_sampled().keys())]),
                              transpose([hstack([s.weight,s.lnl,s.x]) for s in slik.sample()])))))
    else:
        _pool = pool or Pool(nchains)
        try:
            ans = Chains(_pool.map(run_chain,[(main,args,kwargs)]*nchains))
        finally:
            if not pool: _pool.terminate()
        return ans



class SlikPlugin(SlikDict):
    
    def __init__(self,*args,**kwargs):
        """
        Initialization code here.
        """
        super(SlikPlugin,self).__init__(*args,**kwargs)
        
    
    def __call__(self,*args,**kwargs):
        """
        Calculation code here. 
        """
        raise NotImplementedError()


class SlikSampler(SlikDict):
    
    def sample(self,lnl):
        raise NotImplementedError()
    
    
def get_plugin(name):
    """
    Return a :class:`SlikPlugin` class for a given plugin `name`
    
    This::
    
        get_plugin("likelihoods.X")
        
    is equivalent to::
    
        from cosmoslik_plugins.likelihoods.X import X
    
    """
    modname = name.split('.')[-1]
    cls = __import__('cosmoslik_plugins.'+name,fromlist=modname).__getattribute__(modname)
    if not inspect.isclass(cls) or (SlikPlugin not in inspect.getmro(cls) and SlikSampler not in inspect.getmro(cls)):
        raise Exception("Can't load plugin '%s'. It does not appear to be a CosmoSlik plugin."%name)
    return cls
        
        
def get_all_plugins():
    """
    Get all valid CosmoSlik plugins found in the namespace package
    `cosmoslik_plugins`.
    
    Returns:
        dict: mapping of `{class:name}` for each plugin.
    """
    from cosmoslik_plugins import likelihoods, models, samplers, misc
    
    plugins = dict()
    for t in [likelihoods, models, samplers, misc]:
        for _,fullname,_ in  pkgutil.walk_packages(t.__path__,t.__name__+'.',onerror=(lambda x: None)):
            try:
                modname = fullname.split('.')[-1]
                mod = __import__(fullname,fromlist=modname)
                cls = mod.__getattribute__(modname)
                mro = inspect.getmro(cls)
                if SlikPlugin in mro and len(plugins.get(cls,(fullname,))[0])>=len(fullname): 
                    plugins[cls] = fullname
            except Exception:
                pass
        
    return plugins
    

def plugin_getter(module):
    class plugin_getter(object):
        def __init__(self):
            for k,v in get_all_plugins().items():
                vs = v.split(".")
                if vs[1]==module:
                    setattr(self,vs[2],k)
        
        def __getattribute__(self,x):
            try:
                return super().__getattribute__(x)
            except AttributeError: pass
            # try again even though this will fail so we see the import error
            # if it is in fact a real plugin 
            return get_plugin("%s.%s"%(module,x)) 
                
    return plugin_getter()




def lsum(*args):
    """
    If your likelihood function is the sum of a bunch of others, you can do::   

        def __call__(self):
            return lsum(lambda: myfunc1(), lambda: myfunc2())

    This returns the sum ``myfunc1()+myfunc2()``, but never evaluates ``myfunc2()``
    if ``myfunc1()`` returns ``inf``. 

    See also :func:`lsumk` to automatically store the results of ``myfunc1/2``.
    """
    s = 0
    for x in args:
        s+=x()
        if s==inf: break
    return s


def lsumk(lnls,args):
    """
    If your likelihood function is the sum of a bunch of others, you can do::

        def __call__(self):
            self['lnls']={}
            return lsumk(self['lnls'], [('key1',lambda: myfunc1()),
                                        ('key2',lambda: myfunc2())])

    This returns the sum ``myfunc1()+myfunc2()``, but never evaluates
    ``myfunc2()`` if ``myfunc1()`` returns `inf`.  It also stores the result
    ``myfunc1/2`` to ``lnls['key1/2']`` (stores `nan` if a function wasn't
    called)
    """
    s = 0
    for k,f in args:
        if s==inf: lnls[k]=nan
        else:
            lnls[k]=f()
            s+=lnls[k]
    return s



def param_shortcut(*args):
    """
    
    Create a :class:`param` constructor which splices in keyword arguments for you. 

    Example::
    
        param = param_shortcut('start','scale')
        param(1,2)
        > {'start':1, 'scale':2}
    """
    
    class param_shortcut(param):
        def __init__(self,*args2,**kwargs2):
            kwargs2.update(dict(list(zip(args,args2))))
            super(param_shortcut,self).__init__(**kwargs2)
            
    return param_shortcut



def all_kw(ls,exclusions=None):
    for k in (['self','args','kwargs'] 
              + (exclusions if exclusions is not None else [])):
        ls.pop(k,None)
    return ls
    
    
    
def arguments(exclude=None, exclude_self=True, include_kwargs=True, ifset=False):
    """
    Return a dictionary of the current function's arguments (excluding self, and
    optionally other user specified ones or default ones)
    
    Args:
        exclude (list): argument names to exclude
        exclude_self (bool): exclude the "self" argument
        include_kwargs (bool): include the args captured by kwargs
        ifset (bool): exclude keyword arguments whose value is their default
        
    Example::
    
        def foo(x, y=2, **kwargs):
            return arguments()
            
        f(1)        # returns {'x':1, 'y':2}
        f(1,3)      # returns {'x':1, 'y':3}
        f(1,z=5)    # returns {'x':1, 'y':3, 'z':5}
    
    Adapted from: http://kbyanc.blogspot.fr/2007/07/python-aggregating-function-arguments.html
    """
    from inspect import getargvalues, stack, signature
    _, kwname, args = getargvalues(stack()[1][0])[-3:]
    if include_kwargs:
        args.update(args.pop(kwname, []))
    else:
        args.pop(kwname, [])
    if exclude_self: args.pop("self",None)
    if exclude: 
        for e in exclude: args.pop(e,None)
    args.pop('__class__',None) #don't really understand why this shows up sometimes...
    if ifset:
        raise NotImplementedError()
    return args
