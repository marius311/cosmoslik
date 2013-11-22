from cosmoslik import SlikPlugin
from numpy import hstack, arange, array, cumsum, frombuffer, dot, ix_, ones, sqrt, load
from scipy.linalg import cho_factor, cho_solve, inv
import cPickle

class camspec_slik(SlikPlugin):
    
    def __init__(self,
                 lrange=None,
                 like_file=None,
                 dx=None,
                 **kwargs):
        
        super(camspec_slik,self).__init__(**kwargs)
        
        self.labels = [(100,100),(143,143),(217,217),(143,217)]        
        self.freqs = {x:x for x in self.labels}
        self.in_lrange = [(50,1201),(50,2001),(500,2501),(500,2501)]
        self.nl = [u-l for l,u in self.in_lrange]
        self.ells = hstack([arange(*r) for r in self.in_lrange])
        self.out_lrange = self.in_lrange if lrange is None else lrange
        
        for ri,ro in zip(self.in_lrange, self.out_lrange):
            if ro and (ro[0]<ri[0] or ro[1]>ri[1]): 
                raise Exception("Camspec lrange outside of available range.")
        
        self.slice = hstack([arange(*(array(r)+s-l)) 
                             for s,(l,u),r in zip([0]+list(cumsum(self.nl)[:-1]),
                                                  self.in_lrange,
                                                  self.out_lrange)
                             if r])
        
        with open(like_file,'r') as f: self.x, cv = cPickle.load(f)
        if dx is not None: self.x += dx
        todl = self.ells*(self.ells+1)
        self.x *= todl
        cv = ((cv*todl).T*todl).T
        self.cho_cov = inv(cv[ix_(self.slice,self.slice)])
        
        
    def __call__(self, cmb, egfs):
        
        dcl = self.get_x()[self.slice] - hstack(self.get_cl_model(cmb, egfs))[self.slice]
        return dot(dcl,dot(self.cho_cov,dcl))/2
    
    def get_x(self):
#        cal = [p.get(('camspec','cal%i'%i),1) for i in range(3)]
        return self.x #* hstack([ones(n)*a for n,a in zip(self.nl,[cal[0],cal[1],cal[2],sqrt(cal[1]*cal[2])])])
    
    def get_cl_model(self, cmb, egfs):
        return [cmb['cl_TT'][slice(*r)] + 
                egfs(spectra='cl_TT',lmax=r[1],freqs=self.freqs[l])[slice(*r)] 
                for r,l in zip(self.in_lrange, self.labels)]
        
