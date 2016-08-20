from cosmoslik import SlikPlugin, param
from numpy import *
from scipy.linalg import cho_factor, cho_solve, inv
import pickle

class camspec_slik(SlikPlugin):
    
    def __init__(self,
                 lrange=None,
                 like_file=None,
                 dx=None,
                 cal=False,
                 **kwargs):
        
        super(camspec_slik,self).__init__(**kwargs)
        
        self.labels = [(100,100),(143,143),(217,217),(143,217)]        
        self.freqs = {x:x for x in self.labels}
        self.in_lrange = [(50,1201),(50,2001),(500,2501),(500,2501)]
        self.nl = [u-l for l,u in self.in_lrange]
        self.ells = hstack([arange(*r) for r in self.in_lrange])
        self.out_lrange = self.in_lrange if lrange is None else lrange

        if cal:
            self.cal0 = param(start=1,scale=0.0004,range=(0.98,1.02),gaussian_prior=(1.0006,0.0004))
            self.cal2 = param(start=1,scale=0.0015,range=(0.95,1.05),gaussian_prior=(0.9966,0.0015))
        
        for ri,ro in zip(self.in_lrange, self.out_lrange):
            if ro and (ro[0]<ri[0] or ro[1]>ri[1]): 
                raise Exception("Camspec lrange outside of available range.")
        
        self.slice = hstack([arange(*(array(r)+s-l)) 
                             for s,(l,u),r in zip([0]+list(cumsum(self.nl)[:-1]),
                                                  self.in_lrange,
                                                  self.out_lrange)
                             if r])
        
        with open(like_file,'r') as f: self.x, cv = pickle.load(f)
        if dx is not None: self.x += dx
        todl = self.ells*(self.ells+1)
        self.x *= todl
        self.cv = ((cv*todl).T*todl).T
        self.cho_cov = inv(self.cv[ix_(self.slice,self.slice)])
        
        
    def __call__(self, cmb, egfs):
        
        dcl = self.get_x()[self.slice] - hstack(self.get_cl_model(cmb, egfs))[self.slice]
        return dot(dcl,dot(self.cho_cov,dcl))/2
    
    def get_x(self):
        cal = [self.get('cal%i'%i,1) for i in range(3)]
        return self.x * hstack([ones(n)*a for n,a in zip(self.nl,[cal[0],cal[1],cal[2],sqrt(cal[1]*cal[2])])])
    
    def get_cl_model(self, cmb, egfs):
        return [cmb['cl_TT'][slice(*r)] + 
                egfs(spectra='cl_TT',lmax=r[1],freqs=self.freqs[l])[slice(*r)] 
                for r,l in zip(self.in_lrange, self.labels)]
        
    def plot(self, cmb_result, egfs_result, residuals=True, bindl=50, fig=None):

        from matplotlib.pyplot import figure

        if fig is None: fig=figure()
        fig.set_size_inches(6,4*6/1.6)
        
        sl=array(cumsum([0]+self.nl))[:-1]
        ses=array(list(zip(sl,sl+self.nl))) + (array(self.out_lrange) - self.in_lrange)
        def bin(x,auto2d=True):
            bx = array([x[i:i+bindl].mean(axis=0) for i in arange(0,x.shape[0],bindl)[:-1]])
            if auto2d and x.ndim==2: return bin(bx.T,auto2d=False).T
            else: return bx    
        clmodel = hstack(self.get_cl_model(cmb_result,egfs_result))

        for i,((s,e),lbl) in enumerate(zip(ses,['100','143','217','143x217']),1):
            ax=fig.add_subplot(4,1,i)
            
            ax.set_title(lbl)
            if residuals:
                ax.plot(self.ells[s:e],zeros(e-s),c='k')
            else:
                ax.errorbar(self.ells[s:e],clmodel[s:e],c='k')
        
            ax.errorbar(bin(self.ells[s:e]),
                     bin(self.get_x()[s:e]-(clmodel[s:e] if residuals else 0)),
                     yerr=sqrt(diag(bin(self.cv[s:e,s:e]))),ls='',marker='.')
            ax.set_xlim(0,2500)
            if residuals: ax.set_ylim(-100,100)
            else: ax.set_ylim(0,6000)
