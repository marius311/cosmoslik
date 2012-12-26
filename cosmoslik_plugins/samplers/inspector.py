from cosmoslik.plugins import Sampler
from cosmoslik.chains import load_chain, Chains
import cosmoslik 

class inspector(Sampler):    
    """
    
    =============
    The Inspector
    =============
    
    This isn't meant as a true sampler, instead its meant to 
    inspect chains you've already run. You can use it to
    
    - make plots 
    - troubleshoot chains
    - debug modules
    
    Making Plots
    ============
    
    Inspector freezes the chain at a given sample, and returns the 
    internal state of cosmoslik exactly as it was at that point when
    the chain was being run. The internal state is just a
    dictionary, which is returned by::
    
        import cosmoslik
        p = cosmoslik.inspect('params.ini')
    
    After this command, ``p`` contains a dictionary which has all 
    of the keys in the inifile. Sampled parameters have their values
    updated to the given step in the chain. A few extra keys are also
    present. You can access the likelihood modules by::
    
        p['_likelihoods']
        
    Most modules have a ``.plot()`` command. See individual module
    documentation.
     
    """

    def init(self, p):
        c = load_chain(p.pop('output_file'))
        if isinstance(c,Chains): c=c.join()
        self.bestfit = c.best_fit()
        
    def sample(self, x, lnl, p):
        x = [self.bestfit['.'.join(k)] for k in p.get_all_sampled()]
        lnl, p = lnl(x, p)
        yield x, 1, lnl, p


def inspect(params, **kwargs):
    return cosmoslik.sample(params,samplers='inspector', **kwargs).next()[-1]
