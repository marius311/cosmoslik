from cosmoslik_plugins.models.egfs import egfs
from numpy import arange, hstack, loadtxt, zeros, exp, sqrt, pi, mean
import os

class camspec(egfs):
    
    def __init__(self,
                 norm_ell=3000,
                 tsz_norm_fr=143,
                 **kwargs): 
        
        super(camspec,self).__init__(**kwargs)
        
        self.norm_ell = float(norm_ell)
        self.tsz_norm_fr = float(tsz_norm_fr)
        
        self.dir = os.path.dirname(os.path.abspath(__file__))
        padding = zeros(10000)
        def norm(dl,norm=self.norm_ell): return dl/dl[norm]
        self.tsz_template = norm(hstack([[0,0],loadtxt(os.path.join(self.dir,"camspec_templates/tsz_143_eps0.50.dat"))[:,1],padding]))
        self.ksz_template = norm(hstack([[0,0],loadtxt(os.path.join(self.dir,"camspec_templates/cl_ksz_148_trac.dat"))[:,1],padding]))
        self.tszxcib_template = norm(hstack([[0,0],loadtxt(os.path.join(self.dir,"camspec_templates/sz_x_cib_template.dat"))[:,1],padding]))

    def get_colors(self,p):
        return {'ps_100':'orange','ps_143':'orange','ps_217':'orange','ps_143_217':'orange',
                'cib_143':'g','cib_217':'g','cib_143_217':'g',
                'tsz':'m','ksz':'cyan',
                'tsz_cib':'brown'}

    def get_egfs(self, 
                 a_ps_100, a_ps_143, a_ps_217, r_ps, 
                 a_cib_143, a_cib_217, r_cib, n_cib,
                 a_tsz, a_ksz, xi,
                 spectra, lmax, freqs, **kwargs):
        
        freqs = tuple(fr if isinstance(fr,dict) else {'tsz':fr} for fr in freqs)
        
        frlbl = tuple([n for n,(l,u) in {100:(80,120),
                                         143:(130,160), 
                                         217:(200,240)}.items() 
                       if l<mean(fr.values())<u][0] for fr in freqs)
        
        comps={}
        
        comps['tsz'] = a_tsz * self.tsz_template[:lmax] * tszdep(freqs[0]['tsz'],freqs[1]['tsz'], self.tsz_norm_fr) 
        comps['ksz'] = a_ksz * self.ksz_template[:lmax] 
        
        ell = arange(lmax)/self.norm_ell
        
        if frlbl==(100,100):
            comps['ps_100'] = a_ps_100*ell**2
        elif frlbl==(143,143):
            comps['ps_143'] = a_ps_143*ell**2
            comps['cib_143'] = a_cib_143*ell**0.8
            comps['tsz_cib'] = - 2 * xi * sqrt(a_tsz * tszdep(143,143,self.tsz_norm_fr) * a_cib_143) * self.tszxcib_template[:lmax]
        elif frlbl==(217,217):
            comps['ps_217'] = a_ps_217*ell**2
            comps['cib_217'] = a_cib_217*ell**0.8
        elif tuple(sorted(frlbl))==(143,217):
            comps['ps_143_217'] = r_ps*sqrt(a_ps_143*a_ps_217)*ell**2
            comps['cib_143_217'] = r_cib*sqrt(a_cib_143*a_cib_217)*ell**n_cib
            comps['tsz_cib'] = - xi * sqrt(a_tsz * tszdep(143,143,self.tsz_norm_fr) * a_cib_217) * self.tszxcib_template[:lmax]
            
        return comps
            
        
def tszdep(fr1,fr2,fr0):
    """The tSZ frequency dependence."""
    t1,t2,t0 = map(lambda fr: (lambda x0: x0*(exp(x0)+1)/(exp(x0)-1) - 4)(fr/56.78),[fr1,fr2,fr0])
    return t1*t2/t0**2
        
