from cosmoslik.plugins import Model
from numpy import arange, loadtxt, hstack, pi, exp, zeros, ndarray
import os

class egfs(Model):
    """
    
    =======================================================
    How to Create Your Own Extra-Galactic Foreground Models
    =======================================================
       
    To create your own extra-galactic foreground model, create a subclass
    of ``cosmoslik.plugins.models.egfs.egfs`` and override the function ``get_egfs`` to return a 
    dictionary of extra-galactic foreground components. 
    
    Also passed to the ``get_egfs`` function is information from the dataset, such as 
    
    - ``spectra`` : e.g. `cl_TT` or `cl_EE`
    - ``freq`` : a dictionary for different effective frequencies, e.g. 
      ``{'dust': 153, 'radio': 151, 'tsz':150}``
    - ``fluxcut`` : the fluxcut in mJy
    - ``lmax`` : the necessary maximum l
    
    Here's an example egfs model ::
    
        from cosmoslik.plugins.models.egfs import egfs
        
        class MyEgfs(egfs):
        
            def get_egfs(self, p, spectra, freq, fluxcut, lmax, **kwargs):
                return {'single_component': p['amp'] * ones(lmax)}
    
    """
    
    def get_egfs(self, p, *args, **kwargs):
        raise NotImplementedError()
    
    def get_colors(self, p):
        return None
    
    def get(self, p, required):
        
        def get_egfs(*args, **kwargs):
            
            comps = self.get_egfs(p, *args, **kwargs)
            
            if 'plot' in kwargs:
                from matplotlib.pyplot import subplot
                ax = kwargs.pop('ax',None) or subplot(111)
                colors = self.get_colors(p)
                for comp in (lambda key: comps if key is True else key)(kwargs.pop('plot')):
                    ax.plot(comps[comp],label=comp,**({'color':colors[comp]} if colors else {}))
                    
            return sum(comps.values())

        get_egfs.__reduce_ex__ = lambda _: (_dont_pickle,(),None,None,None)
        
        return {'egfs':get_egfs}
    
    
def _dont_pickle(): pass 

