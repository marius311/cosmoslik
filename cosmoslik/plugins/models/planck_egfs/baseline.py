from cosmoslik.plugins import Model
from .. egfs import egfs
from numpy import arange, loadtxt, hstack, pi, exp, zeros, ndarray
import os

class baseline(egfs):
    """
    
    """
    
    def init(self, p): 
        
        self.dir = os.path.dirname(os.path.abspath(__file__))
        padding = zeros(10000)
        def todl(cl,norm=3000): return (lambda dl: dl/dl[norm])((lambda l: cl*l*(l+1)/2/pi)(arange(cl.shape[0])))
        self.clustered_template =  todl(hstack([[0,0],loadtxt(os.path.join(self.dir,"clustered_150.dat"))[:,1],padding]),norm=1500)
        self.tsz_template = todl(hstack([[0,0],loadtxt(os.path.join(self.dir,"tsz.dat"))[:,1],padding]))
        self.ksz_template = todl(hstack([[0,0],loadtxt(os.path.join(self.dir,"ksz_ov.dat"))[:,1],padding]))
        
        if p.get(('egfs','tied_dusty_alpha'),False):
            if ('egfs','dgcl.alpha') in p and ('egfs','dgpo.alpha') in p:
                raise Exception("When setting tied_dusty_alpha=True delete egfs.dgcl.alpha") 

    def get_colors(self, p):
        return {'dgpo':'g','dgcl':'g','radio':'orange','tsz':'magenta','ksz':'cyan'}

    def get_egfs(self, p, spectra, fluxcut, freqs, lmax, **kwargs):
        if spectra != 'cl_TT': return zeros(lmax)
        
        p_egfs = p.get('egfs',{})
        
        if p_egfs.get('tied_dusty_alpha',False): p_egfs['dgcl','alpha'] = p_egfs['dgpo','alpha']
        
        fr1, fr2 = freqs
        
        tilted_clust = (lambda tilt: hstack([self.clustered_template[:1500]*(1500/3000.)**tilt,(arange(1500,20001)/3000.)**tilt]))(p_egfs['dgcl','tilt'])
        
        comps = {'dgpo': p_egfs['dgpo','amp'] * (arange(lmax)/3000.)**2 * plaw_dep(fr1['dust'], fr2['dust'], p_egfs['dgpo','norm_fr'], p_egfs['dgpo','alpha']),
                 'dgcl': p_egfs['dgcl','amp'] * tilted_clust[:lmax] * plaw_dep(fr1['dust'], fr2['dust'], p_egfs['dgcl','norm_fr'], p_egfs['dgcl','alpha']),
                 'radio': p_egfs['radio','amp'] * (fluxcut / p_egfs['radio','norm_fluxcut']) ** (2+p_egfs['radio','gamma']) * (arange(lmax)/3000.)**2 * plaw_dep(fr1['radio'], fr2['radio'], p_egfs['radio','norm_fr'], p_egfs['radio','alpha']),
                 'tsz': p_egfs['tsz','amp'] * self.tsz_template[:lmax] * tszdep(fr1['tsz'],fr2['tsz'],p_egfs['tsz','norm_fr']),
                 'ksz': p_egfs['ksz','amp'] * self.ksz_template[:lmax]}
        
        return comps
    
    
def dBdT(fr1,fr0):
    """ dB/dT at T_CMB """
    dBdT,dBdT0 = map((lambda fr: (lambda x0: x0**4 * exp(x0) / (exp(x0)-1)**2)(fr/57.78)),[fr1,fr0])
    return dBdT/dBdT0  
  
def tszdep(fr1,fr2,fr0):
    """The tSZ frequency dependence."""
    t1,t2,t0 = map(lambda fr: (lambda x0: x0*(exp(x0)+1)/(exp(x0)-1) - 4)(fr/56.78),[fr1,fr2,fr0])
    return t1*t2/t0**2

def plaw_dep(fr1,fr2,fr0,alpha):
    """A power-law frequency dependence."""
    return (fr1*fr2/fr0**2.)**alpha / dBdT(fr1,fr0) / dBdT(fr2,fr0)

def _unpicklable(): pass 

