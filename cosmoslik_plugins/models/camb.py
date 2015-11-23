from __future__ import absolute_import
from cosmoslik import SlikPlugin

name_mapping = {}

class camb(SlikPlugin):
    """
    """
    def __init__(self):
        super(camb,self).__init__()
        import camb as _camb
        self._camb = _camb

    def __call__(self,
                 **params):
        
        params = {name_mapping.get(k,k):v for k,v in params.items()}
        cp = self._camb.set_params(**params)
        result = self._camb.get_results(cp)
        return dict(zip(['cl_%s'%x for x in ['TT','EE','BB','TE']],
                        1e12*result.get_cmb_power_spectra(spectra=['total'])['total'].T))
