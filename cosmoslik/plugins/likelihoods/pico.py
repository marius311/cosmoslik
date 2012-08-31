from cosmoslik.plugins import Likelihood
from numpy import atleast_1d, inf
import pypico

class pico(Likelihood):
    """
    
    ===============
    PICO Likelihood
    ===============
    
    This plugin uses PICO to calculate a likelihood instead of using the actual likelihood.
    
    Parameters
    ---------- 
    
    .. describe:: [pico_lnl].datafile
    
        The PICO datafile to load
        
    .. describe:: [pico_lnl].likelihoods 
    
        The names of the PICO outputs which should be added together to produce the likelihood.
    
    .. describe:: [pico_lnl].on_fail
    
        What to do if we get a ``CantUsePICO`` error (because the parameters are out of bounds).
        Should be one of ``['fail','force','inf']`` to allow CosmoSlik to crash,
        to force PICO to run anyway, or to return infinite log-likelihood, respectively. 
    
    """

    def init(self,p):
        try: datafile = p['pico_lnl','datafile']
        except KeyError: raise Exception("Please specify [pico_lnl]{datafile = ... }")
        else: self.pico = pypico.load_pico(datafile)
        
        self.on_fail = p['pico_lnl'].get('on_fail','fail')
        on_fail_opts = ['fail','force','inf']
        if self.on_fail not in on_fail_opts:
            raise Exception('pico_lnl.on_fail should be one of %s'%on_fail_opts)
            
    def lnl(self,p,model):
        force = (self.on_fail == 'force')
        
        try:
            r = self.pico.get(outputs=atleast_1d(p['pico_lnl']['likelihoods']),force=force,**p)
        except pypico.CantUsePICO:
            if self.on_fail=='inf': r={None:inf}
            else: raise
            
        return sum(r.values())
    
