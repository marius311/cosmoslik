from .egfs import egfs
from numpy import arange

class clust_poisson_egfs(egfs):
    
    def get_colors(self):
        return {'ps':'g','cl':'b'}
        
    def get_egfs(self, 
                 Aps, Acib, ncib,
                 spectra, lmax,  
                 **kwargs):
        
        if spectra != 'cl_TT': return {}
        else:
            return {'ps':  Aps * (arange(lmax)/3000.)**2,
                    'cib': Acib * (arange(lmax)/3000.)**ncib}
