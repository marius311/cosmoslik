from cosmoslik.plugins import Likelihood
from numpy import arange, pi
from cosmoslik.plugins import SubprocessExtension

class xfaster(Likelihood):
    
    def get_required_models(self, p):
        return ['cl_TT','cl_TE','cl_EE','cl_BB','egfs']
    
    def init(self, p):
        self.X = SubprocessExtension('xfaster_likelihood_mod',globals())
        self.X.xfaster_init(p['xfaster','file'])
        self.fluxcut = p['xfaster','fluxcut']
        self.freqs = p['xfaster','freqs']
        
    def lnl(self, p, model):
        
        def to_cl(dl): 
            dl[2:]/=(lambda l: l*(l+1)/2/pi)(arange(2,len(dl)))
            return dl
        
        xf_model = {}
        for x in ['cl_TT','cl_TE','cl_EE','cl_BB']:
            xf_model[x.lower()] = to_cl(model[x] + model['egfs'](x, fluxcut=self.fluxcut, freqs=self.freqs, lmax=len(model[x])))[2:]
        
        return self.X.xfaster_lnl(**xf_model)
