import pkgutil, inspect, sys, os, socket
from multiprocessing import Process, Pipe, Queue
import threading

class CosmoSlikPlugin(object):

    pass


class Likelihood(CosmoSlikPlugin):
    """ A cosmoslik likelihood module. """
    
    def init(self,p):
        """ Do any initialization for the likelihood here. """
        pass

    def lnl(self,p,model):
        """ Return negative log-likelihood. """
        raise NotImplementedError('Implement Likelihood.lnl.')

    def get_extra_params(self,p):
        """ Return a list of nuisance parameter names needed by this likelihood. """
        return {}
        
    def get_required_models(self,p):
        """ Return a list of model components needed by this likelihood, e.g. 'cl_TT' """
        return []
        
        
        
class Model(CosmoSlikPlugin):
    """ A cosmoslik model module. """
    
    def init(self,p):
        """ Do any initialization for the model here. """
        pass

    def get(self,p,required):
        """ Get the model. """
        raise NotImplementedError('Implement Model.get.')
    
    

class Sampler(CosmoSlikPlugin):
    """ A cosmoslik sampler. """
    
    def init(self,p):
        """ Do any initialization for the sampler here. """
        pass
        
    def sample(self, x, lnl, p):
        """ Return a generator which yields samples. """
        raise NotImplementedError('Implement Sampler.sample')
    
    
class Deriver(CosmoSlikPlugin):
    
    def init(self,p):
        """ Do any initialization for the deriver here. """
        
    def add_derived(self,p):
        """ Add derived parameters. """
        raise NotImplementedError('Implement Deriver.add_derived')
    
    
class SubprocessExtension(object):
    """
    This imports a module and runs all of its code in a subprocess.
    Its meant to be used by CosmoSlik plugins which load 
    Python extension modules to facilitate clean exception handling. 
    
    If the plugin loads the extension module via ::
    
        from X import X
        
    then instead use ::
    
        X = SubprocessExtension('X',globals())
    """
    
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
        super_getattr = super(SubprocessExtension,self).__getattribute__
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
    cls = __import__('cosmoslik.plugins.'+name,fromlist=modname).__getattribute__(modname)
    if not inspect.isclass(cls) or CosmoSlikPlugin not in inspect.getmro(cls):
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
    import cosmoslik.plugins
    plugins = dict()
    for _,fullname,_ in  pkgutil.walk_packages(cosmoslik.plugins.__path__,cosmoslik.plugins.__name__+'.'):
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


