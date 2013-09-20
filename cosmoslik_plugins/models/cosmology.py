from cosmoslik import SlikPlugin, all_kw

class cosmology(SlikPlugin):    
    
    def __init__(self,
                 As = 2.4e-9,
                 ns = 0.96,
                 ombh2 = 0.0225,
                 omch2 = 0.12,
                 tau = 0.09,
                 H0 = 70,
                 massive_neutrinos = 3.046,
                 massless_neutrinos = 0.000,
                 omk = 0,
                 omnuh2 = 0,
                 nrun = 0,
                 w = -1,
                 **kwargs):
        
        super(cosmology,self).__init__(**dict(all_kw(locals()),**kwargs))
