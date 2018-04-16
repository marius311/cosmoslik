from cosmoslik import SlikPlugin, arguments
from numpy import hstack, zeros, arange, pi, inf, nan
import os.path as osp

__all__ = ['clik']

# keep a cache of loaded clik files since for some, clik dies if you try to
# reload them. this won't work if you have the same plugin in different folders,
# but its better than nothing. 
loaded_cliks = dict()

# the order in which clik puts the spectra and returns them in get_lmax()
clik_specs = ['TT','EE','BB','TE','TB','EB']

class clik(SlikPlugin):
    """

    """

    def __init__(self,
                 clik_file,
                 auto_reject_errors=False,
                 lensing=None,
                 **nuisance):

        super().__init__(**arguments())

        from clik import clik, clik_lensing
        
        if lensing or "lensing" in clik_file:
            loaded_cliks[clik_file] = self.clik = loaded_cliks[clik_file] if clik_file in loaded_cliks else clik_lensing(clik_file)
            self.lensing = True
            self.clik_specs = ['pp'] + clik_specs
            
        else:
            loaded_cliks[clik_file] = self.clik = loaded_cliks[clik_file] if clik_file in loaded_cliks else clik(clik_file)
            self.lensing = False
            self.clik_specs = clik_specs
                


    def __call__(self, cmb):

        for x, lmax in zip(self.clik_specs,self.clik.get_lmax()):
            if (lmax!=-1) and ((x not in cmb) or (cmb[x].size<lmax+1)):
                raise ValueError("Need the %s spectrum to at least lmax=%i for clik file '%s'."%(x,lmax+1,osp.basename(self.clik_file)))

        cl = hstack([tocl(cmb[x][:lmax+1]) for x, lmax in zip(self.clik_specs,self.clik.get_lmax()) if lmax!=-1])
        nuisance = [self[k] for k in self.clik.get_extra_parameter_names()]
        try:
            lnl = -self.clik(hstack([cl,nuisance]))[0]
            if lnl==0: return inf
            else: return lnl
        except Exception as e:
            if self.auto_reject_errors: return inf
            else: raise

def tocl(dl):
    return hstack([zeros(2),dl[2:]/arange(2,dl.size)/(arange(2,dl.size)+1)*2*pi])
