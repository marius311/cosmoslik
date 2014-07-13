from collections import OrderedDict
import copy, threading, pkgutil, inspect, sys, os, socket
from multiprocessing import Process, Pipe
from numpy import inf, hstack, transpose
import imp, hashlib, time

__all__ = ['load_script','Slik','SlikFunction',
           'SlikDict','SlikPlugin','SlikSampler','param','param_shortcut',
           'SubprocessExtension','get_plugin','get_all_plugins',
           'lsum','all_kw','run_chain']

""" Loaded datafiles will reside in this empty module. """
datafile_module = 'cosmoslik.scripts'
sys.modules[datafile_module]=imp.new_module(datafile_module)



def load_script(scriptfile,*args,**kwargs):   
    """ 
    Read a CosmoSlik script. 
    Passes *args and **kwargs to the script's __init__.
    """ 
    script_module = hashlib.md5(os.path.abspath(scriptfile) + time.ctime()).hexdigest()
    modname='%s.%s'%(datafile_module,script_module)
    mod = imp.load_source(modname,scriptfile)
    sys.modules[modname]=mod
    return Slik(mod.main(*args,**kwargs))


class Slik(object):
    
    def __init__(self,params):
        self.params = params
        
#        def add_slik_functions(slikdict):
#            for k in dir(slikdict): 
#                v=getattr(slikdict,k)
#                if isinstance(v,SlikDict): add_slik_functions(v)
#                elif hasattr(v,'_slik_function'): 
#                    def _make_slik_function(v):
#                        def _slik_function(*args,**kwargs):
#                            return v(self,*args,**kwargs)
#                        _slik_function.__doc__ = v.__doc__
#                        return _slik_function
#                    setattr(self,k,_make_slik_function(v))
#        
#        add_slik_functions(self.params)
        
        
        self._sampled = params.find_sampled()
        
        self.sampler = self.params.sampler
        del self.params.sampler
        
        
    def get_sampled(self):
        return self._sampled
        
        
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
            
        return params(), params
            
            
    def sample(self):
        return self.sampler.sample(self.evaluate)
    
    


        
        
def SlikFunction(func):
    func._slik_function = True
    return func


class SlikDict(dict):
    
    def __init__(self,*args,**kwargs):
        super(SlikDict, self).__init__(*args, **kwargs)
        self.__dict__ = self
        
    def __setitem__(self,k,v):
        if isinstance(k,str):
            setattr(reduce(getattr, k.split('.')[:-1], self), k.split('.')[-1], v)
        else:
            raise ValueError("Parameter key must be string, not %s"%type(k))

    def deepcopy(self):
        cself = copy.copy(self)
        for k,v in vars(self).items():
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
            for k,v in self.__dict__.iteritems(): 
                if isinstance(v,SlikDict): 
                    all_sampled.update(_find_sampled(v,root=root+[k]))
                elif isinstance(v,param):
                    all_sampled[('.'.join(root+[k]))]=v 
            return all_sampled
        
        return OrderedDict(sorted([(k,v) for k,v in _find_sampled(self,[]).iteritems()]))


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
    Runs a CosmoSlik chain, or if nchains!=1 runs a set of chains in parallel (trivially). 

    Non-trivial parallelization, e.g. using MPI or with communication amongst chains, 
    is handled by each cosmoslik.Sampler module. See individual documentation. 


    Arguments:
    ----------
    main - the class object for your top-level SlikPlugin 
    pool - any worker pool which has a pool.map function. 
           default: multiprocessing.Pool(nchains)
    nchains - the number of chains to run
    args - args to pass to main(*args)
    kwargs - kwargs to pass to main(**kwargs)


    Returns:
    --------
    A cosmoslik.chain.Chain instance if nchains=1, 
    otherwise a cosmoslik.chain.Chains instance

    """
    from chains import Chain, Chains
    from multiprocessing.pool import Pool

    if args is None: args=[]
    if kwargs is None: kwargs={}

    if nchains==1:
        slik = Slik(main(*args,**kwargs))
        return Chain(dict(zip(hstack(['weight','lnl',slik.get_sampled().keys()]),
                              transpose([hstack([s.weight,s.lnl,s.x]) for s in slik.sample()]))))
    else:
        _pool = pool or Pool(nchains)
        try:
            ans = Chains(_pool.map(run_chain,[main]*nchains))
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
    
    
no_subproc = False
    
def SubprocessExtension(module_name,globals):
    """
    This imports a module and runs all of its code in a subprocess.
    Its meant to be used by CosmoSlik plugins which load 
    Python extension modules to facilitate clean exception handling. 
    
    If the plugin loads the extension module via ::
    
        from X import X
        
    then instead use ::
    
        X = SubprocessExtension('X',globals())
    """
    if no_subproc:
        exec ('from %s import %s \n'
              'mod = %s')%(module_name,module_name,module_name) in globals, locals()
        return mod
    else:
        return _SubprocessExtension(module_name,globals)
    
    
class _SubprocessExtension(object):
    
    def check_recv(self):
        r = [None]
        def recv(): r[0] = self._conn.recv()
        th = threading.Thread(target=recv)
        th.daemon=True
        th.start()

        while th.is_alive():
            th.join(1e-3)
            if not self._proc.is_alive(): 
                raise Exception('Extension module died.\n%s'%self._subproc_stdout.read())
        
        if isinstance(r[0],Exception): raise r[0]
        else: return r[0]
    
    def _flush_subproc_stdout(self):
        try: self._subproc_stdout.read()
        except IOError: pass
        
    def __init__(self, module_name, globals):
        def subproc_code(conn, fd_out):
            try:
                #redirect output to a pipe back to the main process
                for s in [1,2]: os.dup2(fd_out,s)
                exec ('from %s import %s \n'
                      'mod = %s')%(module_name,module_name,module_name) in globals, locals()
                attrs = {k:(getattr(v,'__doc__',None),callable(v)) for k,v in vars(mod).items()}
                conn.send(attrs)
                while True:
                    (attr, args, kwargs) = conn.recv()
                    if attr in attrs and attrs[attr][1]: conn.send(getattr(mod,attr)(*args, **kwargs))
                    else: conn.send(getattr(mod,attr))
            except Exception, e: 
                conn.send(e)
                raise
                return
    
        self._conn, conn_child = Pipe()
        out_socket_parent, out_socket_child = socket.socketpair()
        out_socket_parent.settimeout(0)
        self._subproc_stdout = out_socket_parent.makefile()
        self._proc = Process(target=subproc_code,args=(conn_child,out_socket_child.fileno()))
        self._proc.daemon = True
        self._proc.start()
        self._attrs = self.check_recv()
        
    def __getattribute__(self, attr):
        super_getattr = super(_SubprocessExtension,self).__getattribute__
        try:
            return super_getattr(attr)
        except AttributeError:
            if attr in self._attrs and self._attrs[attr][1]:
                def wrap_method(*args, **kwargs):
                    self._flush_subproc_stdout()
                    self._conn.send((attr,args,kwargs))
                    return self.check_recv()
                wrap_method.__doc__ = self._attrs[attr][0]
                return wrap_method
            else:
                super_getattr('_conn').send((attr,None,None))
                return super_getattr('check_recv')()

            
def get_plugin(name):
    """
    Return a SlikPlugin class for plugin name. 
    name should be module path relative to cosmoslik.plugins
    """
    modname = name.split('.')[-1]
    cls = __import__('cosmoslik_plugins.'+name,fromlist=modname).__getattribute__(modname)
    if not inspect.isclass(cls) or (SlikPlugin not in inspect.getmro(cls) and SlikSampler not in inspect.getmro(cls)):
        raise Exception("Can't load plugin '%s'. It does not appear to be a CosmoSlik plugin."%name)
    return cls
        
        
def get_all_plugins():
    """
    Gets all valid CosmoSlik plugins found in the 
    namespace package cosmoslik_plugins.
    
    The return value is a dictionary of {class:name} for each plugin.
    
    Valid plugins are any module X which has an attribute X which 
    is a subclass of SlikPlugin. If multiple references to
    X exist in the package, only the shallowest one is returned.  
    """
    import cosmoslik_plugins
    plugins = dict()
    for _,fullname,_ in  pkgutil.walk_packages(cosmoslik_plugins.__path__,cosmoslik_plugins.__name__+'.'):
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


def lsum(*args):
    """
    *args is a list of no argument lambda functions
    returns the sum of the returned values when each is evaluated
    when/if a function returns inf, no other functions are evaluted
    """
    s = 0
    for x in args:
        s+=x()
        if s==inf: break
    return s


def param_shortcut(*args):
    """
    *args is a list of keyword names, e.g. ['start','scale']
    
    returns a param constructor which converts the given *args to **kwargs 
    based on the keywords given above.

    Example:
    
    >> param = param_shortcut('start','scale')
    >> p = param(1,2)
    >> p
    {'start':1, 'scale':2}

    """
    
    class param_shortcut(param):
        def __init__(self,*args2,**kwargs2):
            kwargs2.update(dict(zip(args,args2)))
            super(param_shortcut,self).__init__(**kwargs2)
            
    return param_shortcut



def all_kw(ls,exclusions=None):
    for k in (['self','args','kwargs'] 
              + (exclusions if exclusions is not None else [])):
        ls.pop(k,None)
    return ls

