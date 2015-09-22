from numpy import delete, array, fromstring, loadtxt, dot, arange, diag, hstack, zeros
from scipy.linalg import cho_factor, cho_solve
from cosmoslik import SlikPlugin
import os.path as osp
from itertools import takewhile

class spt_lowl(SlikPlugin):
    
    def __init__(self,
                 which=None,
                 lmin=None,
                 lmax=None,
                 drop=None,
                 **kwargs):
        
        super(spt_lowl,self).__init__(**kwargs)
        
        if which=='s12':
            newdat_file = 'data/s12/spt_lps12_20120828/Spectrum_spt2500deg2_lps12_alternativeCalibrationImplementation.newdat'
        elif which=='k11':
            newdat_file = 'data/k11/bandpowers/Spectrum_spt20082009.newdat'
        else:
            raise ValueError("spt_lowl: 'which' must be one of ['s12','k11']")
        
        newdat_file = osp.join(osp.dirname(__file__),newdat_file)
        
        #Load spectrum and covariance
        with open(newdat_file) as f:
            window_dir = osp.dirname(f.readline())
            while 'TT' not in f.readline(): pass
            self.spec=array([fromstring(f.readline(),sep=' ')[1] for _ in range(47)])
            self.sigma=array([fromstring(f.readline(),sep=' ') for _ in range(94)])[47:]
            
        #Load windows
        self.windows = [loadtxt(osp.join(osp.dirname(newdat_file),'windows',window_dir,'window_%i'%i))[:,1] for i in range(1,48)]
        
        if lmin is not None:
            bmin = sum(1 for _ in takewhile(lambda x: x<lmin, (sum(1 for _ in takewhile(lambda x: abs(x)<.001,w) ) for w in self.windows)))
        else: bmin = 0

        if lmax is not None:
            bmax = sum(1 for _ in takewhile(lambda x: x<lmax, [3251 - sum(1 for _ in takewhile(lambda x: abs(x)<.001,reversed(w)) ) for w in self.windows]))
        else: bmax = 48

        self.spec = self.spec[bmin:bmax]
        self.sigma = self.sigma[bmin:bmax,bmin:bmax]
        self.windows = self.windows[bmin:bmax]
        
        if drop is not None:
        	self.spec = delete(self.spec,drop)
        	self.sigma = delete(delete(self.sigma,drop,0),drop,1)
        	self.windows = delete(self.windows,drop,0)
        
        self.cho_sigma = cho_factor(self.sigma)
        
        self.windowrange = (lambda x: slice(min(x),max(x)+1))(loadtxt(osp.join(osp.dirname(newdat_file),'windows',window_dir,'window_1'))[:,0])
        self.lmax = self.windowrange.stop
        self.ells = array([dot(arange(10000)[self.windowrange],w) for w in self.windows])

        self.freq = {'dust':154, 'radio': 151, 'tsz':153}
        self.fluxcut = 50
        

    

    def __call__(self, 
                 cmb, 
                 egfs):
            
        cl = self.get_cl_model(cmb,egfs)
            
        #Apply windows and calculate likelihood
        dcl = self.get('cal',1)*self.spec-cl
        return dot(dcl,cho_solve(self.cho_sigma, dcl))/2


    def plot(self,p, cl=None,
             ax=None,fig=None,
             residuals=False,
             show_comps=True,
             comps_kw={}):
        
        if cl==None: cl = self.get_cl_model(p['_model'])
        
        if ax is None:
            if fig==None: 
                from matplotlib.pyplot import figure
                fig=figure()
            ax=fig.add_subplot(111)
            
        if not residuals:
            ax.errorbar(self.ells,self.spec,yerr=diag(self.sigma[0]),fmt='.',label='SPT K11')
            ax.plot(self.ells,cl)
            if show_comps: 
                ax.plot(p['_model']['cl_TT'],c='b')
                p['_model']['egfs']('cl_TT', lmax=self.lmax, freqs=(self.freq,self.freq), fluxcut=self.fluxcut, plot=True, ax=ax, **comps_kw)

        else:
            ax.errorbar(self.ells,self.spec-cl,yerr=diag(self.sigma[0]),fmt='.',label='SPT K11')
            ax.plot([self.ells[0],self.ells[-1]],[0]*2)


    def get_cl_model(self, cmb, egfs):
        #Get CMB + foreground model
        cl = (hstack([cmb['cl_TT'],zeros(self.lmax)])[:self.lmax] + 
              egfs(spectra='cl_TT', lmax=self.lmax, freqs=(self.freq,self.freq), fluxcut=self.fluxcut))
        
        #Apply window functions
        return array([dot(cl[self.windowrange],w) for w in self.windows])
    
    




