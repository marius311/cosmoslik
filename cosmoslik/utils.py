import cosmoslik_plugins as _cosmoslik_plugins
import pkgutil as _pkgutil
import inspect as _inspect
import sys as _sys

for _,fullname,_ in  _pkgutil.walk_packages(_cosmoslik_plugins.__path__,_cosmoslik_plugins.__name__+'.'):
    try:
        modname = fullname.split('.')[-1]
        mod = __import__(fullname,fromlist=modname)
        for k,v in vars(mod).items():
            if hasattr(v,'_slik_function'): 
                setattr(_sys.modules[__name__],k,v)
    except Exception as e:
        pass
del modname, mod, fullname, _, k, v
