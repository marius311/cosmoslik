from __future__ import absolute_import
from cosmoslik import SlikPlugin
from numpy import arange, pi


class classy(SlikPlugin):
    """
    Plugin for CLASS.

    Credit: Brent Follin, Teresa Hamill
    """

    #{cosmoslik name : class name}
    name_mapping = {'As':'A_s',
                    'ns':'n_s',
                    'r':'r',
                    'nt':'n_t',
                    'ombh2':'omega_b',
                    'omch2':'omega_cdm',
                    'omnuh2':'omega_ncdm',
                    'tau':'tau_reio',
                    'H0':'H0',
                    'massive_neutrinos':'N_ncdm',
                    'massless_neutrinos':'N_ur',
                    'Yp':'YHe',
                    'pivot_scalar':'k_pivot'}


    def __init__(self):
        super(classy,self).__init__()

        try:
            from classy import Class
        except ImportError:
            raise Exception("Failed to import CLASS python wrapper 'Classy'.")

        self.model = Class()


    def __call__(self,
                 ombh2,
                 omch2,
                 H0,
                 As,
                 ns,
                 tau,
                 omnuh2=0.006,
                 w=None,
                 r=None,
                 nrun=None,
                 omk=0,
                 Yp=None,
                 Tcmb=2.7255,
                 massive_neutrinos=1,
                 massless_neutrinos=2.046,
                 l_max_scalar=3000,
                 l_max_tensor=3000,
                 pivot_scalar=0.002,
                 outputs=[],
                 **kwargs):


        
        self.model.set(output='tCl, lCl, pCl',
                       lensing='yes',
                       l_max_scalars=l_max_scalar,
                       **{self.name_mapping[k]:v for k,v in locals().items() 
                          if k in self.name_mapping and v is not None})
        self.model.compute()

        ell = arange(l_max_scalar+1)
        self.cmb_result = {'cl_%s'%x:(self.model.lensed_cl(l_max_scalar)[x.lower()])*Tcmb**2*1e12*ell*(ell+1)/2/pi
                           for x in ['TT','TE','EE','BB','PP','TP']}
        
        return self.cmb_result

    def get_bao_observables(self, z):
        return {'H':self.model.Hubble(z),
                'D_A':self.model.angular_distance(z),
                'c':1.0,
                'r_d':(self.model.get_current_derived_parameters(['rs_rec']))['rs_rec']}
