
from cosmoslik import SlikPlugin, arguments

class camb(SlikPlugin):
    """
    Compute the CMB power spectrum with CAMB.
    """
    
    #{cosmoslik name : pycamb name}
    name_mapping = {
        'theta':'cosmomc_theta',
        'Yp':'YHe',
        'Neff':'nnu',
    }
    
    def __init__(self,**defaults):
        """
        defaults : dict
            any of the parameters accepted by __call__. these will be their defaults 
            unless explicitly passed to __call__. 
        """
        super().__init__()
        import camb as _camb
        self._camb = _camb
        self.defaults = defaults


    def convert_params(self,**params):
        """
        Convert from CosmoSlik params to pycamb
        """
        params = {self.name_mapping.get(k,k):v for k,v in params.items()}
        if 'cosmomc_theta' in params:
            params['H0'] = None
        return params
        

    def __call__(self,
                 ALens=None,
                 As=None,
                 DoLensing=None,
                 H0=None,
                 k_eta_max_scalar=None,
                 lmax=None,
                 massive_neutrinos=None,
                 massless_neutrinos=None,
                 mnu=None,
                 Neff=None,
                 NonLinear=None,
                 ns=None,
                 ombh2=None,
                 omch2=None,
                 omk=None,
                 pivot_scalar=None,
                 tau=None,
                 theta=None,
                 Yp=None,
                 nowarn=False,
                 **kwargs):
        """
        
        Args
        ----
        nowarn : bool
            don't warn about unrecognized parameters which were passed in
        
        Returns : dict
            dictionary of {'TT':array(), 'TE':array(), ...} giving the CMB Dl's in muK^2
            
        """
        
        if not nowarn and kwargs:
            print('Warning: passing unknown parameters to CAMB: '+str(kwargs)+' (set nowarn=True to turn off this message.)')
        
        
        params = dict(self.defaults, 
                      **{k:v for k,v in arguments(include_kwargs=False, exclude=["nowarn"]).items() 
                         if v is not None})
            
        cp = self._camb.set_params(**self.convert_params(**params))
        self.result = self._camb.get_results(cp)
        return dict(list(zip(['TT','EE','BB','TE'],
                        (cp.TCMB*1e6)**2*self.result.get_cmb_power_spectra(spectra=['total'])['total'].T)))
