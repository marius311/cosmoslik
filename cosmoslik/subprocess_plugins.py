import os, socket
from multiprocessing import Process, Pipe
from threading import Thread
import cPickle as pickle

"""
Provides some utilities for automatically running SlikPlugins in a subprocess.
You might want to do this for example if the plugin calls some C code which
might crash and kill your entire process, but you want to recover gracefully
in such a case (e.g. just reject the step). See documentation for the class
decorator `subprocess_class`.
"""



def subprocess_class(cls=None, 
                     auto_restart=True, 
                     serializer=None, 
                     exclude=None, 
                     force_attr=None):
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
    
        * Only the bound method's code is run in the subprocess. If methods
          return other functions, their code is run in the main process.

            def add(x,y):
                return x+y

            @subprocess_class
            class myclass(object):
                def getadder(self):
                    return add

            #1+2 is calculated in the main subprocess here:
            myclass().getadder()(1,2) 

        * If an attribute is callable, it will by default return a
          SubprocessBoundMethod object which can be called. If you want to
          return the object itself instead, pass the attribute name to
          "force_attr", e.g.:

            class callabledict(dict):
                def __call__(self):
                    #do something 

            @subprocess_class
            class myclass(object):
                mydict = callabledict()

            >>> myclass().mydict
            <SubprocessBoundMethod at 0x7f6747a0e150>


            @subprocess_class(force_attr=['mydict'])
            class myclass(object):
                mydict = callabledict()

            >>> myclass().mydict
            {}

        * Constructors and functions with arbitrary combinations of
          args/kwargs work transparently
        * Arguments and return values must be pickle'able.
        * Setting fields from the main process, e.g. myinst.x = 5, is not
          currently implemented.

    """

    if serializer is None: serializer = PickleSerializer()
    if exclude is None: exclude = []
    exclude += ['_subprocess_class_handler','deepcopy']
    if force_attr is None: force_attr = ['__dict__']


    def _subprocess_class(cls):

        if not issubclass(cls,object):
            raise ValueError("Can only call subprocess_class on new-style classes.")

        # replace __init__, __getattribute__, and __call__ with our own
        # versions which dispatch these to the subprocess. we do so via
        # monkey-patching so we don't change the type of cls which might
        # screw up usage of super() in the original class definition

        super_getattr = cls.__getattribute__
        def override_method(newmethod):
            methodname = newmethod.__name__
            super_method = getattr(cls,methodname)
            def methodwrapper(self,*args,**kwargs):
                if super_getattr(self,'_subprocess_class_handler') is not None:
                    return newmethod(self,*args,**kwargs)
                else:
                    return super_method(self,*args,**kwargs)
            setattr(cls,methodname,methodwrapper)
            return methodwrapper

        @override_method
        def __init__(self,*args,**kwargs):
            self._subprocess_class_handler = SubprocessClassHandler(
                cls,args,kwargs,auto_restart=auto_restart,serializer=serializer,force_attr=force_attr
            )
        @override_method
        def __getattribute__(self,attr):
            if (attr in exclude or 
                not isinstance(self._subprocess_class_handler,SubprocessClassHandler)):
                return super_getattr(self,attr)
            else:
                return self._subprocess_class_handler.get_subprocess_attribute(attr)
        @override_method
        def __call__(self,*args,**kwargs):
            return __getattribute__(self,'__call__')(*args,**kwargs)


        # we've chosen (somewhat arbitrarily) that copy creates a brand new
        # subprocess, so override __copy__ do so
        def deepcopy(self):
            return cls(*self._subprocess_class_handler.args, **self._subprocess_class_handler.kwargs)
        cls.deepcopy = deepcopy


        cls._subprocess_class_handler = True

        return cls

    if cls is not None: return _subprocess_class(cls)
    else: return _subprocess_class


class SubprocessBoundMethod(object):
    def __init__(self,func):
        self.func=func
    def __call__(self,*args,**kwargs):
        self.subc._conn.send((self.func,args,kwargs))
        return self.subc.check_recv()

class SubprocessClassDied(Exception): pass

class SubprocessException(object): 
    def __init__(self,e):
        self.e=e

class SubprocessInitException(SubprocessException): pass

class PickleSerializer(object):
    """
    Serialize using Python's pickle (protocol=1 by default since 2 doesn't
    work with nested SlikPlugin's for some reason)
    """
    def __init__(self,protocol=1):
        self.protocol=protocol
    def load(self,data):
        return pickle.loads(data)
    def dump(self,data):
        return pickle.dumps(data,protocol=self.protocol)

class SubprocessClassHandler(object):
    
    def check_recv(self):
        r = [None]
        def recv(): 
            try: r[0] = self._conn.recv()
            except EOFError: pass
        th = Thread(target=recv)
        th.daemon=True
        th.start()

        while th.is_alive():
            th.join(1e-3)
            if not self._proc.is_alive(): 
                e = SubprocessClassDied('Extension module died.\n%s'%self._subproc_stdout.read())
                if self.auto_restart:
                    self.start()
                raise e
        
        if isinstance(r[0],SubprocessException):
            raise r[0].e

        return r[0]
    
    def _flush_subproc_stdout(self):
        try: self._subproc_stdout.read()
        except IOError: pass
        
    def start(self):
        
        def subproc_code(conn, fd_out):

            #redirect output to a pipe back to the main process
            for s in [1,2]: os.dup2(fd_out,s)

            #instantiate the object
            try:
                self.cls._subprocess_class_handler = None
                obj = self.cls(*self.args,**self.kwargs)
            except Exception as e:
                conn.send(SubprocessInitException(e))
            else:
                conn.send(True)

            #main loop to listen for requests
            while True:
                try:
                    request = conn.recv()
                    if isinstance(request,str):
                        if (request not in self.force_attr and callable(getattr(obj,request))):
                            #getting a function
                            conn.send(SubprocessBoundMethod(request))
                        else:
                            #getting an attribute
                            conn.send(getattr(obj,request))
                    elif isinstance(request,tuple):
                        #calling a function
                        func,args,kwargs = request
                        conn.send(getattr(obj,func)(*args,**kwargs))
                except Exception as e:
                    conn.send(SubprocessException(e))
                    

        serializer = self.serializer
        class wrap_conn(object):
            """
            Wrap a socket connection to use our custom serialization
            """
            def __init__(self,conn):
                self.conn=conn
            def send(self,data):
                return self.conn.send(serializer.dump(data))
            def recv(self):
                return serializer.load(self.conn.recv())

        self._conn, conn_child = map(wrap_conn,Pipe())
        out_socket_parent, out_socket_child = socket.socketpair()
        out_socket_parent.settimeout(0)
        self._subproc_stdout = out_socket_parent.makefile()
        self._proc = Process(target=subproc_code,args=(conn_child,out_socket_child.fileno()))
        self._proc.daemon = True
        self._proc.start()
        assert self.check_recv(), "Failed to instantiate object in an unexpected way"

    
    def __init__(self, cls, args, kwargs, auto_restart=True, serializer=None, force_attr=None):
        self.auto_restart = auto_restart
        self.serializer = serializer
        self.cls = cls
        self.args, self.kwargs = args, kwargs
        self.force_attr = force_attr
        self.start()

    def get_subprocess_attribute(self, attr):
        self._conn.send(attr)
        r = self.check_recv()
        if isinstance(r,SubprocessBoundMethod): 
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
        th = Thread(target=recv)
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
