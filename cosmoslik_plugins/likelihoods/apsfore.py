from cosmoslik.plugins import Likelihood
from numpy import hstack, arange, array, cumsum, frombuffer, dot, ix_, ones, sqrt, loadtxt, mean
from scipy.linalg import cho_factor, cho_solve, inv
import cPickle

class apsfore(Likelihood):
    
    def init(self,p):
        self.labels = [(100,100),(143,143),(217,217),(143,217),'TE','TT']        
        self.freqs = {x:x for x in self.labels}
        self.in_lrange = [(50,1200),(50,2000),(500,2500),(500,2500),(50,2500),(50,2500)]
        self.nl = [u-l for l,u in self.in_lrange]
        self.ells = hstack([arange(*r) for r in self.in_lrange])
#        self.out_lrange = p.get(('camspec','lrange'),self.in_lrange)
        
#        for ri,ro in zip(self.in_lrange, self.out_lrange):
#            if ro and (ro[0]<ri[0] or ro[1]>ri[1]): 
#                raise Exception("Camspec lrange outside of available range.")
#        
#        self.slice = hstack([arange(*(array(r)+s-l)) 
#                             for s,(l,u),r in zip([0]+list(cumsum(self.nl)[:-1]),
#                                                  self.in_lrange,
#                                                  self.out_lrange)
#                             if r])
        
        with open(p['camspec','like_file'],'r') as f: self.x, cv = cPickle.load(f)
        for dx in p.get(('camspec','dx'),[]): self.x += loadtxt(dx)
        self.inv_cov = inv(cv)
        
        
    def get_required_models(self, model):
        return ['cl_TT','cl_TE','cl_EE','egfs']
        
    def lnl(self, p, model):
#        import ipdb; ipdb.set_trace()
        dcl = self.get_x(p) - bin(hstack(self.get_cl_model(p, model)))
        return dot(dcl,dot(self.inv_cov,dcl))/2
    
    def get_x(self,p):
        return self.x
#        cal = [p.get(('camspec','cal%i'%i),1) for i in range(3)]
#        return self.x * hstack([ones(n)*a for n,a in zip(self.nl,[cal[0],cal[1],cal[2],sqrt(cal[1]*cal[2])])])
    
    def get_cl_model(self, p, model):
        return [model['cl_TT'][slice(*r)] + 
                model['egfs']('cl_TT',lmax=r[1],freqs=self.freqs[l])[slice(*r)] 
                for r,l in zip(self.in_lrange, self.labels)[:4]] + \
               [model['cl_TE'][slice(*self.in_lrange[4])],
                model['cl_EE'][slice(*self.in_lrange[5])]]
            
#    def plot(self, p, ax=None):
#        from matplotlib.pyplot import figure
#        if ax is None: ax=figure().add_subplot(111)
#        
#        ax.errorbar(arange(*self.lrange),
#                    self.spec[self.lslice],
#                    yerr=self.errorbars[self.lslice])
#        
#        ax.plot(arange(*self.lrange), self.get_cl_model(p, p['_model']))
#        
        
        
dl=25

def _bin(x):
    return array([mean(x[i:i+dl,...],axis=0) for i in range(0,x.shape[0],dl)])

def bin(x):
    if x.ndim==1: return _bin(x)
    else: return _bin(_bin(x).T).T
    