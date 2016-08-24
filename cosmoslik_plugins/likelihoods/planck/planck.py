from cosmoslik import SlikPlugin, arguments, param_shortcut
from .clik import clik

__all__ = ['planck_2015_highl_TT',
           'planck_2015_highl_TTTEEE',
           'planck_2015_lowl_TT',
           'planck_2015_lowl_TEB']


param = param_shortcut('start','scale')

class planck_2015_highl_TT(clik):
    
    def __init__(
        self,
        clik_file,
        A_cib_217        = param(60,  10,     range=(0,200)),
        A_planck         = param(1,   0.0025, range=(0.9,1.1), gaussian_prior=(1,0.0025)),
        A_sz             = param(5,   3,      range=(0,10)),
        calib_100T       = param(1,   0.001,  range=(0,3),     gaussian_prior=(0.999,0.001)),
        calib_217T       = param(1,   0.002,  range=(0,3),     gaussian_prior=(0.995,0.002)),
        cib_index        = -1.3,   
        gal545_A_100     = param(7,   2,      range=(0,50),    gaussian_prior=(7,2)),
        gal545_A_143     = param(9,   2,      range=(0,50),    gaussian_prior=(9,2)),
        gal545_A_143_217 = param(21,  8.5,    range=(0,100),   gaussian_prior=(21,8.5)),
        gal545_A_217     = param(80,  20,     range=(0,400),   gaussian_prior=(80,20)),
        ksz_norm         = param(2,   3,      range=(0,10)),
        ps_A_100_100     = param(250, 30,     range=(0,4000)),
        ps_A_143_143     = param(45,  10,     range=(0,4000)),
        ps_A_143_217     = param(40,  10,     range=(0,4000)),
        ps_A_217_217     = param(90,  15,     range=(0,4000)),
        xi_sz_cib        = param(0.5, 0.3,    range=(0,1)),
        **kwargs
    ):
        super().__init__(**arguments())



class planck_2015_highl_TTTEEE(planck_2015_highl_TT):
    
    def __init__(
        self,
        A_pol             = 1,
        calib_100P        = 1,
        calib_143P        = 1,
        calib_217P        = 1,
        galf_EE_A_100     = param(0.060, 0.012, range=(0,10), gaussian_prior=(0.060, 0.012)),
        galf_EE_A_100_143 = param(0.050, 0.015, range=(0,10), gaussian_prior=(0.050, 0.015)),
        galf_EE_A_100_217 = param(0.110, 0.033, range=(0,10), gaussian_prior=(0.110, 0.033)),
        galf_EE_A_143     = param(0.10,  0.02,  range=(0,10), gaussian_prior=(0.10,  0.02)),
        galf_EE_A_143_217 = param(0.240, 0.048, range=(0,10), gaussian_prior=(0.240, 0.048)),
        galf_EE_A_217     = param(0.72,  0.14,  range=(0,10), gaussian_prior=(0.72,  0.14)),
        galf_EE_index     = -2.4,
        galf_TE_A_100     = param(0.140, 0.042, range=(0,10), gaussian_prior=(0.140, 0.042)),
        galf_TE_A_100_143 = param(0.120, 0.036, range=(0,10), gaussian_prior=(0.120, 0.036)),
        galf_TE_A_100_217 = param(0.30,  0.09,  range=(0,10), gaussian_prior=(0.30,  0.09)),
        galf_TE_A_143     = param(0.240, 0.072, range=(0,10), gaussian_prior=(0.240, 0.072)),
        galf_TE_A_143_217 = param(0.60,  0.018, range=(0,10), gaussian_prior=(0.60,  0.018)),
        galf_TE_A_217     = param(1.80,  0.54,  range=(0,10), gaussian_prior=(1.80,  0.54)),
        galf_TE_index     = -2.4,
        **kwargs
    ):
        super().__init__(**arguments())
        
        # set beam leakages to zero (if the user didn't pass them in)
        for m in range(5):
            for i,j in ['00','01','02','11','12','22']:
                for xi,xj in ['TE','EE']:
                    k = "bleak_epsilon_{m}_{i}{xi}_{j}{xj}".format(m=m,i=i,j=j,xi=xi,xj=xj)
                    if k not in self: self[k] = 0


class planck_2015_lowl_TT(clik):
    
    def __init__(
        self,
        clik_file,
        A_planck = param(1, 0.0025, range=(0.9,1.1), gaussian_prior=(1,0.0025)),
    ):
        super().__init__(**arguments())


class planck_2015_lowl_TEB(clik):
    
    def __init__(
        self,
        clik_file,
        A_planck = param(1, 0.0025, range=(0.9,1.1), gaussian_prior=(1,0.0025)),
    ):
        super().__init__(**arguments())
