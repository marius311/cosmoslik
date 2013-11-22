import mspec as M
from mspec.utils import pairs
from numpy import dot, arange, product, zeros, load, exp, vstack, hstack, ones, loadtxt, array
from scipy.linalg import cho_solve, cholesky
from cosmoslik import SlikPlugin
from collections import defaultdict
import cPickle

class mspec_lnl(SlikPlugin):
    
    def __init__(self,
                 signal,
                 use=None,
                 cleaning=None,
                 egfs_kwargs=None,
                 ):

        super(mspec_lnl,self).__init__()

        self.signal = signal if isinstance(signal,M.PowerSpectra) else M.load_signal(signal) 
        self.use = use
        self.cleaning = cleaning
        self.egfs_kwargs = egfs_kwargs if egfs_kwargs is not None else defaultdict(lambda: {})        

        self.processed_signal = self.process_signal(self.signal)

        self.signal_matrix_cov = self.processed_signal.get_as_matrix(lrange=self.use).cov
        self.signal_matrix_cho_cov = cholesky(self.signal_matrix_cov), False
        
        
        
    def process_signal(self, sig):
        """Get processed signal, appyling calibration, doing linear combination, etc.."""
              
        return sig if self.cleaning is None else sig.lincombo(self.cleaning)


    def get_cl_model(self, cmb, egfs):
        """ 
        Build an Mspec PowerSpectra object which holds CMB + foreground C_ell's
        for all the required frequencies 
        """
        
        _egfs = egfs if isinstance(egfs,dict) else defaultdict(lambda: egfs)
        
        return M.PowerSpectra({(a,b):cmb['cl_%s%s'%(a[0],b[0])][:lmax] + _egfs[(a,b)](spectra='cl_TT',lmax=lmax,spec=(a,b),**self.egfs_kwargs.get((a,b),{}))
                               for (a,b),(lmin,lmax) in self.use.items()})


    
    def __call__(self,
                 cmb,
                 egfs):
        
        cl_model = self.get_cl_model(cmb, egfs).binned(self.signal.binning)
        cl_model_matrix = cl_model.get_as_matrix(lrange=self.use).spec
        signal_matrix_spec = self.process_signal(self.signal).get_as_matrix(lrange=self.use,get_cov=False).spec
        
        
        #Compute the likelihood  
        dcl = cl_model_matrix - signal_matrix_spec
        return dot(dcl,cho_solve(self.signal_matrix_cho_cov,dcl))/2
