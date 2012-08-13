import pkgutil, inspect, cosmoslik.plugins

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


