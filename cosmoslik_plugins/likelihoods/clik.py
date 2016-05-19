from __future__ import absolute_import
from cosmoslik import SlikPlugin
from numpy import hstack, zeros, arange, pi, inf, nan

class clik(SlikPlugin):
    """

    """

    def __init__(self,
                 clik_file,
                 auto_reject_errors=False,
                 **nuisance):

        super(clik,self).__init__(auto_reject_errors=auto_reject_errors,**nuisance)

        import clik as _clik
        self.clik = _clik.clik(clik_file)


    def __call__(self, cmb):

        cl = hstack(tocl(cmb.get('cl_%s'%x,zeros(lmax+1))[:lmax+1])
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


def tocl(dl):
    return hstack([zeros(2),dl[2:]/arange(2,dl.size)/(arange(2,dl.size)+1)*2*pi])
