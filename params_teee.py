from cosmoslik import param_shortcut, lsum, get_plugin, SlikDict, SlikPlugin, Slik
from numpy import identity, exp, inf, arange
import cPickle
import mspec as M      
        
param = param_shortcut('start','scale')

class main(SlikPlugin):
    
    def __init__(self):
        super(SlikPlugin,self).__init__()
    
        self.cosmo = get_plugin('models.cosmology')(
            logA = param(3.2),
            ns = param(0.96),
            ombh2 = param(0.0221),
            omch2 = param(0.12),
            tau = param(0.09,min=0),
            theta = param(0.010413),
            omnuh2 = 0.000645
        )
        
        self.egfs = SlikDict(
            a_ps_100=param(130,min=0),
            f_pol_100=param(0,.1,range=(0,1)),
            a_ps_143=param(50,min=0),
            a_ps_217=param(50,min=0),
            r_ps=param(1,range=(0,1)),
            a_cib_143=param(10,min=0),
            a_cib_217=param(30,min=0),
            r_cib=param(1,range=(0,1)),
            a_tsz=param(5,min=0),
            a_ksz=param(5,min=0),
            xi=param(1,range=(0,1)),
            n_cib=param(0.8,gaussian_prior=(0.8,0.2),scale=0.1)
        )
                
        self.get_cmb = get_plugin('models.camb')()

        self.bbn = get_plugin('models.bbn_consistency')()
        self.hubble_theta = get_plugin('models.hubble_theta')()
        
        self.lowl = get_plugin('models.pico')(
            datafile='data/pico.tailmonty.plancklike.dat'
        )
        
        
        with open("data/mspec.dat") as f: signal = cPickle.load(f)
        signal.binning = M.get_bin_func('wmap')
        self.mspec = get_plugin('likelihoods.mspec_lnl')(
            signal=signal,
            
            cleaning={('T','100'):{('T','100'):1},
                      ('E','100'):{('E','100'):1}
                      #('T','143'):{('T','143'):1},
                      #('T','217c'):{('T','217'):1,
                      #('T','545'):-0.008}
                      },
                                                       
            use={(('T','100'),('E','100'))  :(50,1201),
                 (('E','100'),('E','100'))  :(50,1201)
                 #(('T','100'),('T','100'))  :(50,1201),
                 #(('T','143'),('T','143'))  :(50,2001),
                 #(('T','217c'),('T','217c')):(500,2501),
                 #(('T','143'),('T','217c')) :(500,2501)
                 },
                                                         
            egfs_kwargs={(('T','100'),('T','100'))  :dict(freqs=(100,100)),
                         (('T','143'),('T','143'))  :dict(freqs=(143,143)),
                         (('T','217c'),('T','217c')):dict(freqs=(217,217)),
                         (('T','143'),('T','217c')) :dict(freqs=(143,217))}
        )
        
        self.get_egfs = get_plugin('models.planck_egfs.camspec')()        
    
        self.priors = get_plugin('likelihoods.priors')(self)
    
        self.sampler = get_plugin('samplers.metropolis_hastings')(
             self,
             num_samples=1000000,
             output_file='test_teee.chain',
             proposal_cov='data/slik.covmat',
             proposal_scale=1,
             output_extra_params=['cosmo.Yp','cosmo.H0']
        )
    
    def __call__(self):
        self.cosmo.As = exp(self.cosmo.logA)*1e-10
        self.cosmo.Yp = self.bbn(**self.cosmo)
        self.cosmo.H0 = self.hubble_theta.theta_to_hubble(**self.cosmo)
        
        #self.cmb_result = self.get_cmb(outputs=['cl_TT','cl_TE','cl_EE'],force=True,**self.cosmo)
        self.cmb_result = self.get_cmb(outputs=['cl_TE','cl_EE'],force=True,**self.cosmo)
        egfs = self.get_egfs(**self.egfs)
        self.egfs_result = {(('T','100'),('T','100'))  :egfs,
                            (('T','143'),('T','143'))  :egfs,
                            (('T','217c'),('T','217c')):egfs,
                            (('T','143'),('T','217c')) :egfs,
                            (('T','100'),('E','100'))  :lambda lmax,**kwargs: self.egfs.a_ps_100 * (arange(lmax)/3000.)**2 * self.egfs.f_pol_100,
                            (('E','100'),('E','100'))  :lambda lmax,**kwargs: self.egfs.a_ps_100 * (arange(lmax)/3000.)**2 * self.egfs.f_pol_100**2}

        
        return lsum(lambda: self.priors(self),
                    lambda: sum(self.lowl(outputs=None,force=True,**self.cosmo).values()),
                    lambda: self.mspec(self.cmb_result,
                                       self.egfs_result))
