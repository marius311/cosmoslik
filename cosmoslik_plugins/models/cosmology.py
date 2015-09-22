from cosmoslik import SlikPlugin, all_kw, param_shortcut, SlikFunction
param = param_shortcut('start','scale')

class cosmology(SlikPlugin):    
    

    def __init__(self,
                 model = '',
                 As = 2.4e-9,
                 ns = 0.96,
                 ombh2 = 0.0225,
                 omch2 = 0.12,
                 tau = 0.09,
                 H0 = None,
                 massive_neutrinos = 1,
                 massless_neutrinos = 2.046,
                 omk = 0,
                 omnuh2 = 0.00064,
                 nrun = 0,
                 w = -1,
                 **kwargs):

        if 'lcdm' in model:
            logA = param(3.2,0.03)
            ns = param(0.96,0.006)
            ombh2 = param(0.0221,0.0002)
            omch2 = param(0.12,0.002)
            tau = param(0.09,0.01,min=0)
            theta = param(0.010413,5e-6)
            omnuh2 = 0.000645
        if 'alens' in model: Alens = param(1,0.1)
        if 'neff' in model: massive_neutrinos = param(3,.2)
        if 'yp' in model: Yp = param(.24,0.1)
        if 'mnu' in model: omnuh2 = param(0,0.001,range=(0,1))
        if 'nrun' in model: nrun = param(0,0.01)
        
        super(cosmology,self).__init__(**dict(all_kw(locals()),**kwargs))

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
    return {prefix+k:v for k,v in latex.items()}

