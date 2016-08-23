from cosmoslik import SlikPlugin, arguments, param_shortcut, SlikFunction
param = param_shortcut('start','scale')

class cosmology(SlikPlugin):    
    

    def __init__(self,
                 model = '',
                #  As = 2.4e-9,
                #  ns = 0.96,
                #  ombh2 = 0.0225,
                #  omch2 = 0.12,
                #  tau = 0.09,
                #  H0 = None,
                #  massive_neutrinos = 1,
                #  massless_neutrinos = 2.046,
                #  omk = 0,
                #  omnuh2 = 0.00064,
                #  nrun = 0,
                #  w = -1,
                 **kwargs):
        
        super().__init__(**arguments(exclude=["model"]))
        model = model.lower()

        if 'lcdm' in model:
            self.As     = param(2.1e-9,  0.03e-9)
            self.ns     = param(0.96,    0.006)
            self.ombh2  = param(0.0221,  0.0002)
            self.omch2  = param(0.12,    0.002)
            self.tau    = param(0.09,    0.01,  min=0)
            self.theta  = param(0.01041, 5e-6)
            self.H0     = None
        if 'alens' in model: 
            self.Alens  = param(1,       0.1)
        if 'neff' in model: 
            self.Neff   = param(3,       0.2)
        if 'yp' in model: 
            self.Yp     = param(0.24,    0.1)
        if 'mnu' in model: 
            self.mnu    = param(0,       0.001, range=(0,1))
        if 'nrun' in model: 
            self.nrun   = param(0,       0.01)
        

@SlikFunction
def cosmo_latex(prefix=''):
    latex = {'H0':r'$H_0$',
             'Yp':r'$Y_p$',
             'logA':r'$\log A$',
             'ns':r'$n_s$',
             'ombh2':r'$\Omega_bh^2$',
             'omch2':r'$\Omega_ch^2$',
             'tau':r'$\tau$',
             'theta':r'$\theta$'}
    return {prefix+k:v for k,v in list(latex.items())}
