from numpy import delete, array, fromstring, loadtxt, dot, arange, diag, hstack, zeros, sqrt
from scipy.linalg import cho_factor, cho_solve
from cosmoslik import SlikPlugin, arguments, param
from cosmoslik_plugins.models.egfs import egfs_specs
import os.path as osp
from itertools import takewhile
import tarfile
from textwrap import dedent

class spt_lowl(SlikPlugin):
    """
    The SPT "low-L" likelihood from Keisler et al. 2011 or Story et al. 2012. 
    
    See `spt_lowl.ipynb` for some examples using this plugin. 
    """
    
    def __init__(self,
                 which='s12',
                 lmin=None,
                 lmax=None,
                 drop=None,
                 egfs='default',
                 cal=None,
                 **kwargs):
        """
        Args:
            which: 'k11' or 's12'
            lmin/lmax (int): restrict the data to this :math:`\ell`-range
            drop (slice, int, or array of ints): remove the bins at these indices
            cal (float): calibration parameter, defined to multiply the data
                power spectrum. (only available if which='s12')
            egfs (callable, None, or 'default'): A callable object which will be
                called to compute the extra-galactic foreground model. It should
                accept keyword arguments `lmax` which specifies the length of
                the array to reuturn, and `egfs_specs` which contains information
                like frequency and fluxcut needed for more sophisticated
                foreground models.  If `egfs` is None, its assumed it will be
                set later (you must do so before computing the likelihood). If
                'default' is specified, the default foreground model will be
                used (see :class:`spt_lowl_egfs`)
        """
        
        super().__init__(**arguments())
        
        if self.egfs == 'default':
            self.egfs = spt_lowl_egfs()
        
        if which=='s12':
            tar_filename = "spt_lps12_20120828.tar.gz"
            if cal is None:
                self.cal = param(start=1, scale=0.026, gaussian_prior=(1,0.026))
        elif which=='k11':
            tar_filename = "bandpowers_spt20082009.tar.gz"
            if not ((cal is None) or cal!=1.):
                raise ValueError(dedent("""
                Can't use a calibration parameter with the 'k11' likelihood, it
                already has the calibration marginalized into the covariance.
                """))
            else:
                self.cal = 1
        else:
            raise ValueError("spt_lowl: 'which' must be one of ['s12','k11']")


        with tarfile.open(osp.join(osp.dirname(__file__),"data",tar_filename)) as tf:
            
            newdat_file = next(f for f in tf.getnames() if "newdat" in f)
            
            #Load spectrum and covariance
            with tf.extractfile(newdat_file) as f:
                while 'TT' not in str(f.readline()): pass
                self.spec=array([fromstring(f.readline(),sep=' ')[1] for _ in range(47)])
                self.sigma=array([fromstring(f.readline(),sep=' ') for _ in range(94)])[47:]
            
            #Load windows
            def load_window(i):
                return loadtxt(tf.extractfile(next(f for f in tf.getnames() if f.endswith("window_%i"%i))))
            self.windows = [load_window(i)[:,1] for i in range(1,48)]
            self.windowrange = (lambda x: slice(int(min(x)),int(max(x)+1)))(load_window(1)[:,0])
        
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
        
        
        self.lmax = self.windowrange.stop
        self.ells = array([dot(arange(10000)[self.windowrange],w) for w in self.windows])

        self.egfs_specs = egfs_specs(
            kind = 'TT',
            freqs = ({'dust':154, 'radio': 151, 'tsz':153},)*2,
            fluxcut = 50
        )
    
    

    def __call__(self, 
                 cmb, 
                 egfs=None,
                 cal=None):
        """
        
        Compute the likelihood. 
        
        Args:
            cmb (dict): A dict which has a key "TT" which holds
                the CMB Dl's in muK^2 
            egfs/cal: Override whatever (if anything) was set in :func:`__init__`.  See
                :func:`__init__` for documentation on these arguments. 
        """
        
        self.update(**{k:v for k,v in arguments().items() if v is not None})
        
        cl = self.get_cl_model(cmb,self.egfs)
            
        #Apply windows and calculate likelihood
        dcl = self.cal*self.spec-cl
        return dot(dcl,cho_solve(self.cho_sigma, dcl))/2


    def plot(self,
             ax=None,
             residuals=False,
             show_comps=False):
        """
        Plot the data (multplied by calibration) and the model. 
        
        Calibration and model taken from previous call to :func:`__call__`.
        
        Args:
            ax (axes): matplotlib axes object
            residual (bool): plot residuals to the model
            show_comps (bool): plot CMB and egfs models separately
        """
        
        if ax is None:
            from matplotlib.pyplot import gca
            ax = gca()
            
        label = 'SPT '+self.which.upper()
        
        if not residuals:
            ax.errorbar(self.ells,self.cal*self.spec,yerr=sqrt(diag(self.sigma)),fmt='.',label=label)
            ax.plot(self.ells,self.get_cl_model(self.cmb,self.egfs))
            if show_comps: 
                ax.plot(self.cmb['TT'][:self.lmax])
                ax.plot(self.egfs(lmax=self.lmax,egfs_specs=self.egfs_specs)[:self.lmax])

        else:
            ax.errorbar(self.ells,self.spec-cl,yerr=sqrt(diag(self.sigma)),fmt='.',label=label)
            ax.plot([self.ells[0],self.ells[-1]],[0]*2)

        ax.set_xlabel(r'$\ell$',size=18)
        ax.set_ylabel(r'$D_\ell\,[\rm \mu K^2]$',size=18)

    def get_cl_model(self, cmb, egfs):
        
        if ('TT' not in cmb):
            raise ValueError('Need the TT spectrum for spt_lowl.')
        if cmb['TT'].size < self.lmax:
            raise ValueError('Need the TT spectrum to at least lmax=%i for spt_lowl.'%self.lmax)
        
        #Get CMB + foreground model
        cl = cmb['TT'][:self.lmax] + egfs(lmax=self.lmax, egfs_specs=self.egfs_specs)[:self.lmax]
        
         #Apply window functions
        return array([dot(cl[self.windowrange],w) for w in self.windows])



class spt_lowl_egfs(SlikPlugin):
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
        self.template_sz = loadtxt(osp.join(rootdir,"dl_sz.txt"))[:,1]
        self.template_cl = loadtxt(osp.join(rootdir,"dl_clustered_point_source.txt"))[:,1]
        self.template_ps = loadtxt(osp.join(rootdir,"dl_poisson_point_source.txt"))[:,1]
        
        
    def __call__(self, lmax, **_):
        return self.Asz*self.template_sz[:lmax] + self.Acl*self.template_cl[:lmax] + self.Aps*self.template_ps[:lmax]
