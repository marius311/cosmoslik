import os, os.path as osp
from cosmoslik import *
from cosmoslik import mpi
from collections import OrderedDict
from numpy.random import uniform
#todo: cleanup handling of these two imports c.f. get_all_plugins
from cosmoslik_plugins.samplers import metropolis_hastings
from cosmoslik_plugins.likelihoods.priors import priors as likelihoods_priors


class priors(metropolis_hastings.metropolis_hastings):
    """
    Sample from the prior.
    """
    
    def __init__(self, params, 
                 output_extra_params=None,
                 **kwargs):
    
        if output_extra_params is None: output_extra_params = []
        super().__init__(params,proposal_update=False,**arguments(exclude=['params','proposal_update']))

        self.output_extra_params = OrderedDict([k if isinstance(k,tuple) else (k,dtype('float').name) for k in output_extra_params])
        self.sampled = params.find_sampled()

        self.uniform_priors = {}
        _uniform_priors = likelihoods_priors(params).uniform_priors
        for k in self.sampled:
            pr = [(n,l,u) for (n,l,u) in _uniform_priors if n==k]
            if len(pr)>2:
                raise Exception("Can't sample prior for parameters with multiple priors at once.")
            elif len(pr)==0:
                raise Exception("Prior for parameter '%s' is improper."%k)
            else:
                self.uniform_priors[k] = pr[0][1:]

        # todo: gaussian priors
        
    def draw_x(self):
        return [uniform(*self.uniform_priors[k]) for k in self.sampled]
                
    def _mcmc(self, _, lnl):
    
        for _ in range(int(self.num_samples/float(max(1,mpi.get_size()-1)))):
            cur_x = self.draw_x()
            cur_lnl, cur_extra = lnl(*cur_x)
            yield(sample(cur_x, cur_lnl, 1, cur_extra))
