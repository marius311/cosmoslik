from cosmoslik.plugins import Likelihood
from numpy import load, zeros, dot, arange, sqrt, diag, loadtxt, array
from scipy.linalg import cho_factor, cho_solve

class onespec(Likelihood):
    
    def init(self,p):
        
        p = p['onespec']
        ells, spec = (load(p['spectrum']) if p['spectrum'].endswith('npy') else loadtxt(p['spectrum'])).T
        self.spec = zeros(max(ells)+1)
        self.spec[array(ells,dtype=int)] = spec
        cov = load(p['covariance']) if p['covariance'].endswith('npy') else loadtxt(p['covariance'])
        self.errorbars = sqrt(diag(cov))
        self.lrange = p['lrange']
        self.lslice = slice(*self.lrange)
        self.cho_cov = cho_factor(cov[self.lslice,self.lslice])
        self.egfs_kwargs = p['egfs_kwargs']
        
    def get_required_models(self, model):
        return ['cl_TT','egfs']
        
    def lnl(self, p, model):
        dcl = self.spec[self.lslice] - self.get_cl_model(p, model)
        return dot(dcl,cho_solve(self.cho_cov,dcl))/2
    
    def get_cl_model(self, p, model):
        return model['cl_TT'][self.lslice] + model['egfs']('cl_TT',lmax=self.lrange[1],**self.egfs_kwargs)[self.lslice]
        
    def plot(self, p, ax=None):
        from matplotlib.pyplot import figure
        if ax is None: ax=figure().add_subplot(111)
        
        ax.errorbar(arange(*self.lrange),
                    self.spec[self.lslice],
                    yerr=self.errorbars[self.lslice])
        
        ax.plot(arange(*self.lrange), self.get_cl_model(p, p['_model']))
        
        