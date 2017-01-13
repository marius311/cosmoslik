from cosmoslik import SlikPlugin, param
from numpy import inf

class priors(SlikPlugin):
    
    def __init__(self,params):
        self.gaussian_priors = []
        self.uniform_priors = []

        for k in params.find_sampled(): 
            if hasattr(params[k], "gaussian_prior"): self.add_gaussian_prior(k, *params[k].gaussian_prior)
            if hasattr(params[k], "uniform_prior"): self.add_uniform_prior(k, *params[k].uniform_prior)
            if hasattr(params[k], "range"): self.add_uniform_prior(k, *params[k].range)
            if hasattr(params[k], "min"): self.add_uniform_prior(k, params[k].min,inf)
            if hasattr(params[k], "max"): self.add_uniform_prior(k, -inf, params[k].max)
    
    def add_gaussian_prior(self,param,mean,std):
        self.gaussian_priors.append((param,mean,std))
        
    def add_uniform_prior(self,param,min,max):
        self.uniform_priors.append((param,min,max))
    
    def __call__(self,params):
        param_ok = lambda n: n in params and not isinstance(params[n],param)
        for n, lower, upper in self.uniform_priors:
            if param_ok(n) and not (lower<=params[n]<=upper): 
                return inf
        
        return sum((params[n]-c)**2./2/w**2 for n,c,w in self.gaussian_priors if param_ok(n))
        
