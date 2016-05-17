import threading, os, socket
from multiprocessing import Process, Pipe

"""
Provides some utilities for automatically running SlikPlugins in a subprocess.
You might want to do this for example if the plugin calls some C code which
might crash and kill your entire process, but you want to recover gracefully
in such a case (e.g. just reject the step). See documentation for the class
decorator `subprocess_class`.
"""



def subprocess_class(auto_restart=True):
    """
    A class decorator which makes all of the class's code run in a
    subprocess. You can also specify whether to automatically restart the
    subprocess and reinstantiate the class if the subprocess ever dies.

    Example:

        @subprocess_class(auto_restart=True)
        class myclass(object):
            x=3
            def f(self):
                return 11
            
            def die(self):
                sys.exit()
        

        >>> myinst = myclass()
        >>> myinst.x 
        3
        >>> myinst.f()
        11
        >>> myinst.die()
        Traceback (most recent call last):
        ...
        SubprocessDiedException
        >>> myinst.x  #this works because the subprocess was automatically restarted
        3

    Notes:
    
        * Constructors and functions with arbitrary combinations of
          args/kwargs work transparently
        * Arguments and return values must be pickle'able.
        * Setting fields, e.g. myinst.x = 5 not currently implemented.

    """
    def _subprocess_class(classobj):
        def create_subprocess_class(*args,**kwargs):
            return SubprocessClass(classobj,auto_restart=auto_restart,*args,**kwargs)
        return create_subprocess_class
    return _subprocess_class


class SubprocessInstanceMethod(object):
    def __init__(self,func):
        self.func=func
    def __call__(self,*args,**kwargs):
        self.subc._conn.send((self.func,args,kwargs))
        return self.subc.check_recv()

class SubprocessClassDied(Exception): pass

class SubprocessException(object): 
    def __init__(self,e):
        self.e=e

class SubprocessClass(object):
    
    def check_recv(self):
        r = [None]
        def recv(): r[0] = self._conn.recv()
        th = threading.Thread(target=recv)
        th.daemon=True
        th.start()

        while th.is_alive():
            th.join(1e-3)
            if not self._proc.is_alive(): 
                e = SubprocessClassDied('Extension module died.\n%s'%self._subproc_stdout.read())
                if self.auto_restart:
                    self._start()
                raise e
        
        return r[0]
    
    def _flush_subproc_stdout(self):
        try: self._subproc_stdout.read()
        except IOError: pass
        
    def _start(self):
        
        classobj,args,kwargs = self._start_args
        
        def subproc_code(conn, fd_out):
            #redirect output to a pipe back to the main process
            for s in [1,2]: os.dup2(fd_out,s)

            #instantiate the object
            obj=classobj(*args,**kwargs)

            #listen for requests
            while True:
                try:
                    request = conn.recv()
                    if isinstance(request,str):
                        if callable(getattr(obj,request)):
                            conn.send(SubprocessInstanceMethod(request))
                        else:
                            conn.send(getattr(obj,request))
                    elif isinstance(request,tuple):
                        func,_args,_kwargs = request
                        conn.send(getattr(obj,func)(*_args,**_kwargs))
                except Exception as e:
                    conn.send(SubprocessException(e))
                        
        #launch the subprocesses
        self._conn, conn_child = Pipe()
        out_socket_parent, out_socket_child = socket.socketpair()
        out_socket_parent.settimeout(0)
        self._subproc_stdout = out_socket_parent.makefile()
        self._proc = Process(target=subproc_code,args=(conn_child,out_socket_child.fileno()))
        self._proc.daemon = True
        self._proc.start()
    
    def __init__(self, classobj, auto_restart=True, *args, **kwargs):
        self.auto_restart=auto_restart
        self._start_args = (classobj,args,kwargs)
        self._start()

    def __call__(self,*args,**kwargs):
        return self.__getattribute__('__call__',force_subproc=True)(*args,**kwargs)
        
    def __getattribute__(self, attr, force_subproc=False):
        try:
            if not force_subproc:
                return super(SubprocessClass,self).__getattribute__(attr)
        except AttributeError:
            pass
        
        self._conn.send(attr)
        r = self.check_recv()
        if isinstance(r,SubprocessException):
            raise r.e
        if isinstance(r,SubprocessInstanceMethod): 
            r.subc=self
        return r






# ------------------------------------------------------
# Older code which only wrapped Python extension modules
# ------------------------------------------------------


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
