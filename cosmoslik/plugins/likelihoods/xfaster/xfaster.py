from cosmoslik.plugins import Likelihood
from xfaster_likelihood_mod import xfaster_likelihood_mod as X
from numpy import arange, pi

class xfaster(Likelihood):
    
    def get_required_models(self, p):
        return ['cl_TT','cl_TE','cl_EE','cl_BB','egfs']
    
    def init(self, p):
        X.xfaster_init(p['xfaster','file'])
        self.fluxcut = p['xfaster','fluxcut']
        self.freqs = p['xfaster','freqs']
        
    def lnl(self, p, model):
        
        def to_cl(dl): return dl/(lambda l: l*(l+1)/2/pi)(arange(len(dl)))
        
        xf_model = {}
        for x in ['cl_TT','cl_TE','cl_EE','cl_BB']:
            xf_model[x.lower()] = to_cl(model[x] + model['egfs'](x, fluxcut=self.fluxcut, freqs=self.freqs, lmax=len(model[x])))[2:]
        
        return X.xfaster_lnl(**xf_model)
