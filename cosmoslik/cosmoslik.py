from collections import OrderedDict
import copy, threading, pkgutil, inspect, sys, os, socket
from multiprocessing import Process, Pipe

__all__ = ['load_script','create_script','load_init_script','Slik','SlikFunction',
           'SlikDict','SlikPlugin','SlikSampler','param',
           'SubprocessExtension','get_plugin','get_all_plugins']


def load_script(scriptfile):   
    """ 
    Read a CosmoSlik script. 
    """ 
    import imp
    mymod = imp.new_module('_test')
    with open(scriptfile) as f: code=f.read()
    exec code in mymod.__dict__
    
    #Clean out modules 
    from types import ModuleType
    for k,v in mymod.__dict__.items(): 
        if isinstance(v,ModuleType): mymod.__dict__.pop(k)
        
    #Clean out builtins
    mymod.__dict__.pop('__builtins__')
    
    return Slik(**mymod.__dict__)


def create_script(**kwargs):
    return Slik(**kwargs)

def load_init_script(scriptfile):
    p = load_script(scriptfile)
    p.init_plugins()
    return p



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
        self.sampler.init(self)
        
        
#    def init_plugins(self):
#        for k in ['likelihoods','models','derivers','sampler']:
#            for m in atleast_1d(getattr(self.params,k,[])): 
#                if mpi.is_master(): 
#                    print 'Initializing %s...'%m.__class__.__name__
#                m.init(self.params)
#
#        self.sampler = self.params.sampler
#        del self.params.sampler
#        

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
        
        if len(args)!=0:
            for k,v in zip(self.get_sampled().keys(),args): params[k]=v
        else:
            for k,v in kwargs.items(): params[k]=v
            
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
    
    def get(self,k,default):
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
    
    
no_subproc = True
    
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
    Return a CosmoSlikPlugin class for plugin name. 
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
    namespace package cosmoslik.plugins.
    
    The return value is a list of (name, class, type) for each plugin.
    
    Valid plugins are any module X which has an attribute X which 
    is a subclass of CosmoSlikPlugin. If multiple references to
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
            if CosmoSlikPlugin in mro and len(plugins.get(cls,(fullname,))[0])>=len(fullname): 
                for t in [Likelihood, Model, Sampler, Deriver]: 
                    if t in mro: 
                        typ = t.__name__
                        break
                else: continue
                plugins[cls] = (fullname,typ)
        except Exception:
            pass
        
    plugins = sorted([(fullname,cls,typ) for cls, (fullname, typ) in plugins.items()])
    return plugins



