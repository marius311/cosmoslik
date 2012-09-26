from cosmoslik.plugins import Likelihood
from numpy import *
from itertools import islice, chain
from scipy.linalg import cho_factor, cho_solve
import re

class newdat(Likelihood):
    """
    
    
    """
    
    xs = ['TT','EE','BB','EB','TE','TB']
    
    def init(self, p):
        
        self.binrange = {x:slice(*p.get(('newdat','binrange'),{}).get(x,[None])) for x in self.xs}
        self.fluxcut = p['newdat','fluxcut']
        self.freqs = p['newdat','freqs']
        
        self.egfs = p.get(('newdat','egfs'),'egfs')
        
        fac = p.get(('newdat','factor'),1)
        
        with open(p['newdat','file']) as f:
            name = remove_comments(f.next())
            nxs = map(int,remove_comments(f.next()).split())
            has_calib_uncert, calib, calib_err = remove_comments(f.next()).split()
            has_beam_uncertain, beam, beam_err = remove_comments(f.next()).split()
            ilike = int(remove_comments(f.next()))
            self.bands = {}
            
            self.binslice = hstack([arange(cnx,cnx+nx)[self.binrange[x]] for x, nx, cnx in zip(self.xs,nxs,cumsum(hstack([0,nxs])))])

            for x,nx in zip(self.xs,nxs):
                if nx!=0:
                    if f.next().strip()!=x: raise Exception('Error reading newdat file. Expected bandpowers in order %s'%self.xs)
                    self.bands[x] = array([fromstring(remove_comments(s),sep=' ') for s in islice(f,nx)])
                    self.bands[x][:,1:5] *= fac
                    for _ in islice(f,nx): pass #ignore correlation matrix
        
            self.lmax = max(chain(*[b[:,6] for b in self.bands.values()]))

            self.cov = cho_factor(array([fromstring(remove_comments(s),sep=' ') for s in islice(f,sum(nxs))])[ix_(self.binslice,self.binslice)]*fac**2)            
                
                
    def get_required_models(self, p):
        return ['cl_%s'%x for x in self.bands.keys()]
            
    
    def plot(self,
             ax=None,
             cl=None, 
             p=None, 
             #show_comps=False,
             show_data=True,
             show_model=True,
             show_comps=False,
             residuals=False,
             data_color='k',
             model_color='k'):
        
        if cl is None: cl=self.get_cl_model(p, p['_model'])
        if ax is None: 
            from matplotlib.pyplot import figure
            ax=figure().add_subplot(111)
            
        cl_model = self.get_cl_model(p)
            
        ells = [mean(lrange) for lrange in self.bands['TT'][self.binrange['TT'],[5,6]]]
        if residuals:
            if show_data: ax.errorbar(ells,self.bands['TT'][self.binrange['TT'],1] - cl_model['TT'][self.binrange['TT']],
                                      yerr=self.bands['TT'][self.binrange['TT'],[3,2]].T,
                                      fmt='.',
                                      color=data_color)
        else:
            if show_data: 
                ax.errorbar(ells,self.bands['TT'][self.binrange['TT'],1],
                            yerr=self.bands['TT'][self.binrange['TT'],[3,2]].T,
                            color=data_color,
                            fmt='.')
            if show_model: 
                ax.plot(ells,cl_model['TT'][self.binrange['TT']],color=model_color)
                
            if show_comps:
                ax.plot(p['_model']['cl_TT'],label='CMB')
                p['_model'][self.egfs]('cl_TT',
                                       fluxcut=min(self.fluxcut,self.fluxcut),
                                       freqs=self.freqs,
                                       lmax=self.lmax,
                                       plot=True,
                                       ax=ax)
            
        return ax

    def get_cl_model(self, p, model=None):
        if model is None: model=p['_model']
        
        cl_model = {}
        for x in self.xs:
            if x in self.bands.keys():
                cl_model[x] = model['cl_%s'%x]
                if len(cl_model[x]) < self.lmax: raise Exception('Newdat likelihood needs cl_%s to lmax=%i'%(x,self.lmax))
                cl_model[x] += model[self.egfs]('cl_%s'%x, fluxcut=self.fluxcut, freqs=self.freqs, lmax=len(cl_model[x]))
                
                
        cl_model = {x:array([mean(cl_model[x][lmin:lmax+1]) for (lmin,lmax) in self.bands[x][:,[5,6]]]) for x in self.xs if x in self.bands}

        return cl_model


    def lnl(self, p, model):
        
        cl_model = self.get_cl_model(p, model)
        dcl = hstack([cl_model[x]-self.bands[x][:,1] for x in self.xs if x in self.bands])[self.binslice]
        
        return dot(dcl,cho_solve(self.cov,dcl))/2
        
        
def remove_comments(s):
    return re.sub('#.*','',s).strip()



#  !File Format:
#  !name
#  !n_bands_TT n_EE, n_BB, n_EB, n_TE, n_TB
#  !has_calib_uncertain calib(amplitude) calib_err(power)
#  !has_beam_uncertain beam beam_err
#  !ilike (0: Gaussian, 1: all x-factor, 2: specified have x-factor)
#  !loop over {
#  ! band-types
#  ! band info: num obs + - x_factor l_min l_max use_x
#  ! correlation-matrix (ignored)
#  ! }  
#  ! covariance matrix
            