
from cosmoslik import SlikPlugin, arguments
from numpy import arange, pi


class classy(SlikPlugin):
    """
    Compute the CMB power spectrum with CLASS.

    Based on work by: Brent Follin, Teresa Hamill
    """

    #{cosmoslik name : class name}
    name_mapping = {
        'As':'A_s',
        'lmax':'l_max_scalars',
        'mnu':'m_ncdm',
        'Neff':'N_ncdm',
        'ns':'n_s',
        'nt':'n_t',
        'ombh2':'omega_b',
        'omch2':'omega_cdm',
        'omk':'Omega_k',
        'pivot_scalar':'k_pivot',
        'r':'r',
        'tau':'tau_reio',
        'Tcmb':'T_cmb',
        'Yp':'YHe',
    }


    def __init__(self,**defaults):
        super().__init__()
        from classy import Class
        self.model = Class()
        self.defaults = defaults


    def convert_params(self,**params):
        """
        Convert from CosmoSlik params to CLASS
        """
        params = {self.name_mapping.get(k,k):v for k,v in params.items()}
        if 'theta' in params:
            params['100*theta_s'] = 100*params.pop('theta') 
        params['lensing'] = 'yes' if params.pop('DoLensing',True) else 'no'
        return params
        
        
    def __call__(self,
                 As=None,
                 DoLensing=True,
                 H0=None,
                 lmax=None,
                 mnu=None,
                 Neff=None,
                 nrun=None,
                 ns=None,
                 ombh2=None, 
                 omch2=None,
                 omk=None,
                 output='tCl, lCl, pCl',
                 pivot_scalar=None,
                 r=None,
                 tau=None,
                 Tcmb=2.7255,
                 theta=None,
                 w=None,
                 Yp=None,
                 nowarn=False,
                 **kwargs):
        
        if not nowarn and kwargs:
            print('Warning: passing unknown parameters to CLASS: '+str(kwargs)+' (set nowarn=True to turn off this message.)')
        
        params = dict(self.defaults,**{k:v for k,v in arguments(include_kwargs=False, exclude=["nowarn"]).items() if v is not None})
        self.model.set(self.convert_params(**params))
        self.model.compute()

        lmax = params['lmax']
        ell = arange(lmax+1)
        self.cmb_result = {x:(self.model.lensed_cl(lmax)[x.lower()])*Tcmb**2*1e12*ell*(ell+1)/2/pi
                           for x in ['TT','TE','EE','BB','PP','TP']}

        self.model.struct_cleanup()
        self.model.empty()
        
        return self.cmb_result

    def get_bao_observables(self, z):
        return {'H':self.model.Hubble(z),
                'D_A':self.model.angular_distance(z),
                'c':1.0,
                'r_d':(self.model.get_current_derived_parameters(['rs_rec']))['rs_rec']}
