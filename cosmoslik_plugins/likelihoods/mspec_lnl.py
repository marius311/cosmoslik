from numpy import dot, arange, product, zeros, load, exp, vstack, hstack, ones, loadtxt, array, dot, sqrt, diag
from scipy.linalg import cho_solve, cholesky, inv
from cosmoslik import SlikPlugin, all_kw
from collections import defaultdict
import pickle
from scipy.stats import gamma
from scipy.special import erfinv

class mspec_lnl(SlikPlugin):
    
    def __init__(self,
                 signal,
                 use=None,
                 resgal=None,
                 cleaning=None,
                 egfs_kwargs=None,
                 cal=None,
                 ):

        super(mspec_lnl,self).__init__(**all_kw(locals()))
        
        import mspec as M
        from mspec.utils import pairs


        if not isinstance(signal,M.PowerSpectra): self.signal = M.load_signal(signal) 
        if egfs_kwargs is None: self.egfs_kwargs = defaultdict(lambda: {})    

        processed_signal = self.process_signal(self.signal,docal=False)
        self.signal_matrix_cov = processed_signal.get_as_matrix(lrange=self.use).cov
        self.signal_matrix_cho_cov = cholesky(self.signal_matrix_cov), False
        
        
        
    def process_signal(self, sig, docal=True, docov=True):
        """Get processed signal, appyling calibration, doing linear combination, etc.."""
        
        if self.cleaning: 
            sig = sig.lincombo(self.cleaning,docov=docov)

        if self.resgal: 
            sig = sig - self.resgal

        if docal and self.cal: 
            sig = sig.lincombo({k:{k:self.cal.get('_'.join(k),1)} for k in sig.get_maps()},normalize=False,docov=docov)
            
        return sig


    def get_cl_model(self, cmb, egfs, use=None):
        """ 
        Build an Mspec PowerSpectra object which holds CMB + foreground C_ell's
        for all the required frequencies 
        """
        
        _egfs = egfs if isinstance(egfs,dict) else defaultdict(lambda: egfs)
        
        return M.PowerSpectra({(a,b):cmb['cl_%s%s'%(a[0],b[0])][:lmax] + _egfs[(a,b)](spectra='cl_TT',lmax=lmax,spec=(a,b),**self.egfs_kwargs.get((a,b),{}))
                               for (a,b),(lmin,lmax) in list((use or self.use).items())})


    
    def __call__(self,
                 cmb,
                 egfs):
        
        cl_model = self.get_cl_model(cmb, egfs).binned(self.signal.binning)
        cl_model_matrix = cl_model.get_as_matrix(lrange=self.use).spec
        signal_matrix_spec = self.process_signal(self.signal,docov=False).get_as_matrix(lrange=self.use).spec
        
        #Compute the likelihood  
        dcl = cl_model_matrix - signal_matrix_spec
        return dot(dcl,cho_solve(self.signal_matrix_cho_cov,dcl))/2
        
        
        
        
    """
    Plotting / analyzing tools:
    """
        
        
    def plot(self,spec,cmb,egfs):
        """
        Plot a single residual spectrum
        """
        from matplotlib.pyplot import errorbar, ylim
        bin=self.signal.binning
        s = bin(slice(*self.use[spec]))
        cl_model = self.get_cl_model(cmb, egfs).binned(bin)[spec][s]
        processed_signal_spec = self.process_signal(self.signal,docal=True,docov=False)
        processed_signal_cov  = self.process_signal(self.signal,docal=False,docov=True)
        errorbar(processed_signal_spec.ells[s],
                 processed_signal_spec[spec][s] - cl_model,
                 yerr=sqrt(diag(processed_signal_cov.cov[spec,spec][s,s])),
                 ls='',marker='.')
        ylim(*{'TT':(-50,50),'TE':(-10,10),'EE':(-5,5)}[spec[0][0]+spec[1][0]])        
        
    def plot_all(self,cmb,egfs,which=None,ncol=3,size=4,aspect=1.6):
        """
        Plot many residual spectra
        """
        from matplotlib.pyplot import figure, subplot, plot, xlim, title
        if which is None: use=self.use
        else: use={k:self.use[k] for k in which}

        nrow = len(use)/ncol+1
        fig=figure()
        if size is not None: fig.set_size_inches(size*ncol,size*nrow/aspect)
        fig.subplots_adjust(hspace=0.4)
        
        lmin,lmax=min([v[0] for v in list(use.values())]),max([v[1] for v in list(use.values())])
        
        for i,k in enumerate(sorted(use),1):
            subplot(nrow,ncol,i)
            plot([lmin,lmax],[0]*2,'k')
            self.plot(k,cmb,egfs)
            xlim(lmin,lmax)
            title(k)        


    def chi2(self,cmb,egfs,use=None,nparams=0):
        """
        Return a tuple of (chi2, dof, pte, nsig).
        
        Parameters:
        -----------
        cmb/egfs: The cmb and egfs result.
        use: Which spectra and lranges to include in the chi2.
        nparams: Number to subtract from the d.o.f. for calculating the PTE. 
        """
        if use is None: use=self.use
        use = {k:(lambda x: x if x is not None else self.use[k])(use.get(k)) for k in use}
        cl_model = self.get_cl_model(cmb, egfs).binned(self.signal.binning)
        cl_model_matrix = cl_model.get_as_matrix(lrange=use).spec
        signal_matrix = self.process_signal(self.signal).get_as_matrix(lrange=use)
            
        dcl = cl_model_matrix - signal_matrix.spec
        chi2 = dot(dcl,dot(inv(signal_matrix.cov),dcl))
        k = cl_model_matrix.size - nparams
        pte = gamma.cdf(k/2.,chi2/2.)
        nsig = sqrt(2.)*erfinv(1-pte)
        return (chi2,k,pte,nsig)
