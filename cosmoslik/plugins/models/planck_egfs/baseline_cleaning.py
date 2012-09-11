from cosmoslik.plugins import Model
from .. egfs import egfs
from numpy import arange, loadtxt, hstack, pi, exp, zeros, sqrt
import os

class baseline_cleaning(egfs):
    """
    
    """
    
    def init(self, p): 
        
        self.dir = os.path.dirname(os.path.abspath(__file__))
        padding = zeros(10000)
        self.norm_ell = float(p.get('egfs',{}).get('norm_ell',3000))
        def todl(cl,norm=self.norm_ell): return (lambda dl: dl/dl[norm])((lambda l: cl*l*(l+1)/2/pi)(arange(cl.shape[0])))
        self.clustered_template =  todl(hstack([[0,0],loadtxt(os.path.join(self.dir,"clustered_150.dat"))[:,1],padding]))
        self.tsz_template = todl(hstack([[0,0],loadtxt(os.path.join(self.dir,"tsz.dat"))[:,1],padding]))
        self.ksz_template = todl(hstack([[0,0],loadtxt(os.path.join(self.dir,"ksz_ov.dat"))[:,1],padding]))
        
        if p.get(('egfs','tied_dusty_alpha'),False):
            if ('egfs','dgcl.alpha') in p and ('egfs','dgpo.alpha') in p:
                raise Exception("When setting tied_dusty_alpha=True delete egfs.dgcl.alpha") 

    def get_colors(self, p):
        return {'dgpo':'g','dgcl':'g','radio':'orange','tsz':'magenta','ksz':'cyan'}

    def get_egfs(self, p, spectra, fluxcut, freqs, lmax, **kwargs):
        if spectra != 'cl_TT': return zeros(lmax)
        
        lowp = p.get(('egfs','low_fr'),{})
        highp = p.get(('egfs','high_fr'),{})
        
        if lowp.get('tied_dusty_alpha',False): lowp['dgcl','alpha'] = lowp['dgpo','alpha']
        
        fr1, fr2 = freqs
        
        dustcomp = {}
        
        for i,fr in enumerate(freqs):
            if fr['dust']<300:
                dustcomp[i] = {'dgpo': sqrt(lowp['dgpo','amp']) * (arange(lmax)/3000.) * plaw_dep(fr['dust'], lowp['dgpo','norm_fr'], lowp['dgpo','alpha']),
                               'dgcl_lin': sqrt(lowp['dgcl','amp_lin'] * self.clustered_template[:lmax]) * plaw_dep(fr['dust'], lowp['dgcl','norm_fr'], lowp['dgcl','alpha']),
                               'dgcl_nonlin': sqrt(lowp['dgcl','amp_nonlin'] * (arange(lmax)/self.norm_ell)**p.get(('dgcl','tilt' ),0.8)) * plaw_dep(fr['dust'], lowp['dgcl','norm_fr'], lowp['dgcl','alpha'])}
            else:
                dustcomp[i] = {'dgpo': sqrt(highp['dgpo','amp']) * (arange(lmax)/3000.), 
                               'dgcl_lin': sqrt(highp['dgcl','amp_lin'] * self.clustered_template[:lmax]),
                               'dgcl_nonlin': sqrt(highp['dgcl','amp_nonlin'] * (arange(lmax)/self.norm_ell)**p.get(('dgcl','tilt' ),0.8))}


        ffr1, ffr2 = fr1['dust'], fr2['dust']

        if 130<ffr1<150 and ffr2>300 or ffr1>300 and 130<ffr2<150: corr = highp['dgpo','cor143']
        if 210<ffr1<230 and ffr2>300 or ffr1>300 and 210<ffr2<230: corr = highp['dgpo','cor217']
        else: corr = 1
        
        comps = {}
        for x in ['dgpo','dgcl_lin','dgcl_nonlin']:
            comps[x] = corr*dustcomp[0][x]*dustcomp[1][x]
            
        comps.update({'radio': lowp['radio','amp'] * (fluxcut / lowp['radio','norm_fluxcut']) ** (2+lowp['radio','gamma']) * (arange(lmax)/3000.) * plaw_dep2(fr1['radio'], fr2['radio'], lowp['radio','norm_fr'], lowp['radio','alpha']),
                      'tsz': lowp['tsz','amp'] * self.tsz_template[:lmax] * tszdep(fr1['tsz'],fr2['tsz'],lowp['tsz','norm_fr']),
                      'ksz': lowp['ksz','amp'] * self.ksz_template[:lmax]})
            
        return comps
    
    
def dBdT(fr1,fr0):
    """ dB/dT at T_CMB """
    dBdT,dBdT0 = map((lambda fr: (lambda x0: x0**4 * exp(x0) / (exp(x0)-1)**2)(fr/57.78)),[fr1,fr0])
    return dBdT/dBdT0  
  
def tszdep(fr1,fr2,fr0):
    """The tSZ frequency dependence."""
    t1,t2,t0 = map(lambda fr: (lambda x0: x0*(exp(x0)+1)/(exp(x0)-1) - 4)(fr/56.78),[fr1,fr2,fr0])
    return t1*t2/t0**2

def plaw_dep(fr,fr0,alpha):
    """A power-law frequency dependence."""
    return (fr/fr0)**alpha / dBdT(fr,fr0)

def plaw_dep2(fr1,fr2,fr0,alpha):
    """A power-law frequency dependence."""
    return (fr1*fr2/fr0**2)**alpha / dBdT(fr1,fr0) / dBdT(fr1,fr0)


def _unpicklable(): pass 

