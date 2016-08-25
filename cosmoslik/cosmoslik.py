from collections import OrderedDict
from functools import reduce
from importlib import import_module
from inspect import getmro, getargvalues, stack, getargspec
from numpy import inf, nan, hstack, transpose, isinf
from pkgutil import walk_packages
import argparse
import copy, pkgutil, sys, os
import imp, hashlib, time
from hashlib import md5


__all__ = ['load_script','Slik','SlikFunction',
           'SlikDict','SlikPlugin','SlikSampler','param','param_shortcut','get_all_plugins',
           'lsum','lsumk','all_kw','run_chain','SlikMain', 'arguments', 'get_caller']           

""" Loaded datafiles will reside in this empty module. """
datafile_module = 'cosmoslik.scripts'
sys.modules[datafile_module]=imp.new_module(datafile_module)


def load_script(script):   
    """ 
    Read a CosmoSlik script from a file and return the :class:`SlikPlugin`
    object.
    """
    if isinstance(script,str):
        script_module = md5(open(script).read().encode()).hexdigest()
        modname='%s.%s'%(datafile_module,script_module)
        if modname is sys.modules:
            mod = sys.modules[modname]
        else:
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
    
    argspec = getargspec(main.__init__)
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


class Slik(object):
    
    def __init__(self,params):
        self.params = params
        self._sampled = params.find_sampled()
        self.sampler = self.params.sampler
        del self.params.sampler
        if "priors" not in params:
            from . import likelihoods
            params.priors = likelihoods.priors(params)
        
        
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
            for k,v in kwargs.items(): params[k]=v
        elif len(args)!=len(self.get_sampled()):
            raise ValueError("Expected %i parameters, only got %i."%(len(self.get_sampled()),len(args)))
        else:
            for k,v in zip(self.get_sampled().keys(),args): params[k]=v
            
        if isinf(params.priors(params)):
            # if we've hit any hard priors don't even evaluate this parameter set
            return inf, params
        else:
            # call priors() again since potentially some could depend on derived
            # parameters which are added only by the __call__
            return params()+params.priors(params), params
            
            
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
            *k1, k2 = k.split('.')
            setattr(reduce(getattr, k1, self), k2, v)
        else:
            raise ValueError("Parameter key must be string, not %s"%type(k))

    def deepcopy(self):
        cself = copy.copy(self)
        for k,v in vars(self).items():
            if isinstance(v,SlikDict): 
                setattr(cself,k,v.deepcopy())
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
    
    
    
def get_all_plugins(ignore_import_errors=True):
    """
    Search recursively through the namespace package `cosmoslik_plugins` for any `SlikPlugin`s.
    
    Returns:
        dict: keys are the shortest qualified name for the object and values are a type object. 
    """
    from cosmoslik_plugins import likelihoods, models, samplers, misc
    
    plugins = dict()
    for t in [likelihoods, models, samplers, misc]:
        for _,fullname,_ in  walk_packages(t.__path__,t.__name__+'.',onerror=(lambda x: None)):
            try:
                mod = import_module(fullname)
            except Exception as e:
                if not ignore_import_errors: raise
                plugins[e] = fullname
            else:
                for name in dir(mod):
                    attr = getattr(mod,name)
                    if (attr is not SlikPlugin
                        and attr is not SlikSampler
                        and hasattr(attr,'__mro__')
                        and (SlikPlugin in getattr(attr,'__mro__') 
                             or SlikSampler in getattr(attr,'__mro__'))
                        and (attr not in plugins
                             or len(plugins[attr]) > len(fullname+'.'+name))):
                        
                        plugins[attr] = fullname+'.'+name
        
    return {v:k for k,v in plugins.items()}
    
    


def plugin_getter(module):
    
    plugins = {tuple(k.split('.')):v for k,v in get_all_plugins().items()}
    if isinstance(module,str): module = tuple(module.split('.'))
    
    class _plugin_getter(object):
        
        def __init__(self):
            for k,v in plugins.items():                    
                if k[:len(module)]==module and k!=module:
                    nxt = k[len(module)]
                    setattr(self,nxt,plugin_getter(module+(nxt,)))
        
        def __call__(self,*args,**kwargs):
            return plugins[module+(module[-1],)](*args,**kwargs)
        
        def __getattribute__(self,name):
            val = super().__getattribute__(name)
            if isinstance(val,Exception): raise AttributeError(name) from val
            else: return val

    if module in plugins:
        return plugins[module]
    else:
        return _plugin_getter()

    
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
    
    
    
def get_caller(depth=1):
    """
    Get the calling function. Works in many typical (but not all) cases.
    
    Thanks to: http://stackoverflow.com/a/39079070/1078529
    
    Example::
    
        def foo():
            return get_caller()
            
        foo() #returns the 'foo' function object

    or
            
        def foo():
            return bar()
            
        def bar():
            return get_caller(depth=2) 
            
        foo() #returns the 'foo' function object
    """
    fr = sys._getframe(depth)   # inspect.stack()[1][0]
    co = fr.f_code
    for get in (
        lambda:fr.f_globals[co.co_name],
        lambda:getattr(fr.f_locals['self'], co.co_name),
        lambda:getattr(fr.f_locals['cls'], co.co_name),
        lambda:fr.f_back.f_locals[co.co_name], # nested
        lambda:fr.f_back.f_locals['func'],  # decorators
        lambda:fr.f_back.f_locals['meth'],
        lambda:fr.f_back.f_locals['f'],
        ):
        try:
            func = get()
        except (KeyError, AttributeError):
            pass
        else:
            if func.__code__ == co:
                return func
    raise AttributeError("func not found")
    
    
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
