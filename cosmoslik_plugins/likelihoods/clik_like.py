from cosmoslik.plugins import Likelihood
from numpy import hstack, zeros, arange, pi, inf, nan

class clik_like(Likelihood):
    """
    
    [clik]{
    
        files = {'highL':'path...'}
    
        [highl]{
            nuisance...
        }
    
    }
    
    TODO: implement nuisance parameters
    
    """
    
    def __init__(self,
                 clik_file):
        
        super(clik_like,self).__init__()
        
        import clik
        self.clik = clik.clik(clik_file)
        
        
    def __call__(self, cmb):
        lnl = -self.clik(hstack(tocl(cmb('cl_%s'%x,zeros(lmax+1))[:lmax+1])
                         for x, lmax in zip(['TT','EE','BB','TE','TB','EB'],self.clik.get_lmax()) 
                         if lmax!=-1))[0]
        if lnl==0: return inf
        else: return lnl
            
    
def tocl(dl): 
    return hstack([zeros(2),dl[2:]/arange(2,dl.size)/(arange(2,dl.size)+1)*2*pi])
    