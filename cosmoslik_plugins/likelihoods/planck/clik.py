from cosmoslik import SlikPlugin, arguments
from numpy import hstack, zeros, arange, pi, inf, nan

__all__ = ['clik']

# keep a cache of loaded clik files since for some, clik dies if you try to
# reload them. this won't work if you have the same plugin in different folders,
# but its better than nothing. 
loaded_cliks = dict()

class clik(SlikPlugin):
    """

    """

    def __init__(self,
                 clik_file,
                 auto_reject_errors=False,
                 **nuisance):

        super().__init__(auto_reject_errors=auto_reject_errors,**nuisance)

        from clik import clik
        
        if clik_file in loaded_cliks:
            self.clik = loaded_cliks[clik_file]
        else:
            loaded_cliks[clik_file] = self.clik = clik(clik_file)


    def __call__(self, cmb):

        cl = hstack(tocl(cmb.get(x,zeros(lmax+1)),lmax,x)
                    for x, lmax in zip(['TT','EE','BB','TE','TB','EB'],self.clik.get_lmax())
                    if lmax!=-1)
        nuisance = [self[k] for k in self.clik.get_extra_parameter_names()]
        try:
            lnl = -self.clik(hstack([cl,nuisance]))[0]
            if lnl==0: return inf
            else: return lnl
        except Exception as e:
            if self.auto_reject_errors: return inf
            else: raise

def tocl(dl,lmax,x):
    if lmax is not None:
        if dl.size<lmax+1: raise ValueError("Need the %s spectrum to at least lmax=%i for clik."%(x,lmax+1))
        dl=dl[:lmax+1]
    return hstack([zeros(2),dl[2:]/arange(2,dl.size)/(arange(2,dl.size)+1)*2*pi])
