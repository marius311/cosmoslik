from cosmoslik import SlikPlugin, SlikDict
from numpy import arange, loadtxt, hstack, pi, exp, zeros, ndarray
import os

class egfs_specs(SlikDict):
    """
    This class stores information needed to calculate the extra-galactic foreground
    contribution to some particular power spectrum. This information is, 
    
    kind : 'TT', 'TE', ...
    freqs : tuple of dicts
        the pair of frequencies being correlated. each entry isnt a number,
        rathers its a dict  with keys 'dust', 'radio', and 'tsz', specfying the
        band center for each type of component in GHz
    fluxcut : the fluxcut 
        the fluxcut in mJy
    """


class egfs(SlikPlugin):
    """
    An 
       
    To create your own extra-galactic foreground model, create a subclass
    of ``cosmoslik.plugins.models.egfs.egfs`` and override the function ``get_egfs`` to return a 
    dictionary of extra-galactic foreground components. 
    
    Also passed to the `get_egfs` function is information from the dataset, such as 
    
    - `spectra` : e.g. `cl_TT` or `cl_EE`
    - `freq` : a dictionary for different effective frequencies, e.g. 
      `{'dust': 153, 'radio': 151, 'tsz':150}`
    - `fluxcut` : the fluxcut in mJy
    - `lmax` : the necessary maximum l
    
    Here's an example egfs model:
    
        from cosmoslik.plugins.models.egfs import egfs
        
        class MyEgfs(egfs):
        
            def get_egfs(self, p, spectra, freq, fluxcut, lmax, **kwargs):
                return {'single_component': p['amp'] * ones(lmax)}
    
    """
    
    def get_egfs(self, *args, **kwargs):
        raise NotImplementedError()
    
    def get_colors(self):
        return {}
    
    def plot(self, ax=None, *args, **kwargs):
        from matplotlib.pyplot import gca
        if ax is None: ax = gca()
        comps = self.get_egfs(*args, **kwargs)
        colors = self.get_colors()
        for comp in (lambda key: comps if key is True else key)(kwargs.pop('plot')):
            ax.plot(comps[comp],label=comp,**({'color':colors[comp]} if comp in colors else {}))
    
    def __call__(self, **kwargs):
        
        def get(**kwargs2):
            kwargs2.update(kwargs)
            return sum(self.get_egfs(**kwargs2).values())
        
        return get
            
            
                
