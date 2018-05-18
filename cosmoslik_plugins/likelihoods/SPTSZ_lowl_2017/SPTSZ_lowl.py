from numpy import outer, delete, array,ones, fromstring, loadtxt, dot, arange, diag, hstack, zeros, log, shape
from scipy.linalg import cho_factor, cho_solve
from cosmoslik import SlikPlugin, arguments, param
import os.path as osp
from itertools import takewhile
from numpy.linalg import slogdet
from math import sqrt
import os

class SPTSZ_lowl(SlikPlugin):
    
    def __init__(self,
                 
                 lmin=None,
                 lmax=None,
                 ab_on=False,
                 cal = 1,
                 **kwargs):
        
        self.cal=cal
        self.egfs = SPTSZ_lowl_egfs()
        self.ab_on=ab_on
        newdat_file = 'data/SPTSZ_bandpowers_and_errors/spt_lps12_20120828/Spectrum_spt2500deg2_lps12_alternativeCalibrationImplementation_no_beam_cov.newdat'   
        
        newdat_file = osp.join(osp.dirname(__file__),newdat_file)
            #Load spectrum and covariance
        with open(newdat_file) as f:
            window_dir = osp.dirname(f.readline())
            while 'TT' not in f.readline(): pass
            self.spec=array([fromstring(f.readline(),sep=' ')[1] for _ in range(47)])
            self.sigma=array([fromstring(f.readline(),sep=' ') for _ in range(47)])
        
        self.windows = array([loadtxt(osp.join(osp.dirname(newdat_file),'windows',window_dir,'window_%i'%i))[:,1] for i in range(1,48)])

        self.windowrange = (lambda x: slice(int(min(x)),int(max(x)+1)))(loadtxt(osp.join(osp.dirname(newdat_file),'windows',window_dir,'window_1'))[:,0])
        self.lmax = self.windowrange.stop
        self.ells = array([dot(arange(10000)[self.windowrange],w) for w in self.windows])
        
        if lmin is not None:
            bmin = sum(1 for _ in takewhile(lambda x: x<lmin, (sum(1 for _ in takewhile(lambda x: abs(x)<.001,w) ) for w in self.windows)))
        else: bmin = 0

        if lmax is not None:
            bmax = sum(1 for _ in takewhile(lambda x: x<lmax, [self.lmax - sum(1 for _ in takewhile(lambda x: abs(x)<.001,reversed(w)) ) for w in self.windows]))
        else: bmax = 48
        

        self.beam_corr = self.spt_beam_cor(freqs=[150,150])[self.windowrange,self.windowrange]
        
        if ((lmin is not None) or (lmax is not None)):
            self.spec = self.spec[bmin:bmax]
            self.sigma = self.sigma[bmin:bmax,bmin:bmax]
            self.windows = self.windows[bmin:bmax]
         
    

    def __call__(self, 
                 cmb):
            
        self.db,beam_cov = self.get_cl_model_and_beam_cov(cmb,self.egfs)
        #Apply windows and calculate likelihood
        cho_sigma=cho_factor(beam_cov+self.sigma)
        self.deldb = self.spec -self.db*self.cal
        return dot(self.deldb,cho_solve(cho_sigma, self.deldb))/2. 


    
    def get_cl_model_and_beam_cov(self, cmb, egfs):
        #(Get CMB + foreground model ) * (1+beams_errs)
        self.egfs=egfs(lmax=self.lmax)[:self.lmax]
        self.cmb = cmb['TT'][:self.lmax]
        self.dl = (self.cmb[self.windowrange] + self.egfs[self.windowrange])
        if self.ab_on:
            self.dl*=(1+.26*1.23e-3*self.DlnClDlnl(self.dl))
        #Apply window functions
        return (dot(self.windows,self.dl), dot(self.windows,dot(self.beam_corr*outer(self.dl,self.dl),self.windows.T )))    

    
    def DlnClDlnl(self,y):
        
        x = arange(self.lmax)[self.windowrange]
        lnx=log(x)
        lny=log(y)
        return array([0]+[(lny[i+1]-lny[i])/(lnx[i+1]-lnx[i]) for i in arange(len(y)-1)])
    




    def spt_beam_cor(self,freqs=[150,150],lmax=4000): 


        dir = os.path.dirname(os.path.abspath(__file__)) #abosolute path to current file, keep beam templates in subdirectory
        
        #put beam templates in dictionary labled according to name and year
        names = ['dc','alpha','wobble','xtalk','outer','inner','venus'] 
        freq = '150'
        years = ['2008','2009','2010','2011']
        year_weights = {'2008': 0.065919966,'2009': 0.22589127,'2010': 0.28829848,'2011': 0.41989028}
        beam_errs = {}
        for name in names:
            for year in years:
                    beam_err_file = os.path.join(dir,
                        'data/SPTSZ_beam_errors/errgrid_'+name+'_'+year+'_'+freq+'.txt') 
                    beam_errs[name+year]=loadtxt(beam_err_file)[:,1]
        mode_names = names[:2]+[names[2:][i]+years[j] for i in range(5) for j in range(4)]      
     
        #we only care about the four year average for dc and alpha so lets do that now
        for year in years:
                for k in mode_names[:2]:
                    if year == '2008':
                        beam_errs[k] = beam_errs[k+year]*year_weights[year]                    
                    else:
                        beam_errs[k] += beam_errs[k+year]*year_weights[year]
                        
                 
                 
            
        correlated = zeros([lmax,lmax])
            
        uncorrelated_root_weight = zeros([lmax,lmax])
            
        uncorrelated = zeros([lmax,lmax])       
            
            
            #now we convert the beam errors into dCl/Cl and sum them together according to correlation
     
        for k in mode_names:         
            if k in mode_names[:2]:
                b_err = (1+beam_errs[k][:lmax])**(-2)-1
                correlated += outer(b_err,b_err)
            elif k in mode_names[2:14]:
                b_err = ((1+ beam_errs[k][:lmax])**(-2)-1)*year_weights[k[-4:]]**.5
                uncorrelated_root_weight += outer(b_err,b_err)
            else:
                b_err = ((1+beam_errs[k][:lmax])**(-2)-1)*year_weights[k[-4:]]
                uncorrelated += outer(b_err,b_err)
                        
        
        return correlated+uncorrelated_root_weight+uncorrelated
                      
class SPTSZ_lowl_egfs(SlikPlugin):
    """
    The default foreground model for the SPT low-L likelihood from Story/Keisler
    et al. which features (tSZ+kSZ), CIB clustering, and CIB Poisson templates,
    with free amplitudes and priors on these amplitudes.
    """
    
    def __init__(self,
                 Asz = param(start=5,  scale=3, min=0, gaussian_prior=(5.5,  3.0)),
                 Acl = param(start=20, scale=3, min=0, gaussian_prior=(19.3, 3.5)),
                 Aps = param(start=5,  scale=3, min=0, gaussian_prior=(5,    2.5)),
                 **kwargs):
        super().__init__(**arguments())
        
        rootdir = osp.join(osp.dirname(__file__),"data")
        self.template_sz = loadtxt(osp.join(rootdir,"SPTSZ_lowl_foreground_templates/SZ_template.txt"))[:,1]
        self.template_cl = loadtxt(osp.join(rootdir,"SPTSZ_lowl_foreground_templates/cluster_template.txt"))[:,1]
        self.template_ps = loadtxt(osp.join(rootdir,"SPTSZ_lowl_foreground_templates/poisson_template.txt"))[:,1]
        
        
    def __call__(self, lmax, **_):
        return self.Asz*self.template_sz[:lmax] + self.Acl*self.template_cl[:lmax] + self.Aps*self.template_ps[:lmax]
                    
