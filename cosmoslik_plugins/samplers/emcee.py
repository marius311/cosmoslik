from cosmoslik import *
from cosmoslik import mpi
from cosmoslik.chains.chains import combine_covs
from emcee import EnsembleSampler
from collections import OrderedDict, namedtuple
from numpy import ix_, zeros, ceil
from numpy.random import multivariate_normal

mcmcsample = namedtuple("mcmcsample",["weight","lnl","x"])

class emcee(SlikSampler):
    
    def __init__(
        self,
        params,
        num_samples,
        nwalkers=100,
        proposal_cov=None,
        output_extra_params=None
    ):
    
        if output_extra_params is None: output_extra_params = []
        super().__init__(**arguments(exclude=['params']))

        self.output_extra_params = OrderedDict([k if isinstance(k,tuple) else (k,dtype('float').name) for k in output_extra_params])
        self.sampled = params.find_sampled()
        self.x0 = [params[k].start for k in self.sampled]
        self.proposal_cov = self.initialize_covariance(self.sampled)
        if mpi.get_size()>1:
            self.pool = mpi.get_pool()
        else:
            self.pool = None
    
    def initialize_covariance(self, sampled):
        """
        Prepare the proposal covariance based on anything passed to
        self.proposal_cov, defaulting to the `scale` of each sampled parameter
        otherwise.
        """
        in_covs = [{k:v.scale for k,v in list(sampled.items()) if hasattr(v,'scale')}]
        if self.proposal_cov is not None:
            in_covs += (self.proposal_cov if isinstance(self.proposal_cov,list) else [self.proposal_cov])
        names, covs = combine_covs(*in_covs)
        
        missing = [s for s in sampled if s not in names]
        if missing:
            raise ValueError("Parameters %s not in covariance and no scale given."%missing)
        
        idxs = [names.index(s) for s in sampled]
        return covs[ix_(idxs,idxs)]
                 
    def sample(self, lnl):
        
        # distribute walkers initially according to covariance estimate
        pos = multivariate_normal(
            [v.start for v in self.sampled.values()],
            self.proposal_cov,
            size=self.nwalkers
        )
        
        # diferent sign convention with emcee
        def lnprob(x):
            l, p = lnl(*x)
            return -l, p 
                        
        # step each walker once, yield them one-by-one, repeat
        weight = zeros(self.nwalkers)
        sampler = EnsembleSampler(self.nwalkers, len(self.sampled), lnprob, pool=self.pool)

        for i in range(int(ceil(self.num_samples/self.nwalkers))):
            posnext, prob, state, blobs = sampler.run_mcmc(pos,1)
            for i, (x, xnext, l) in enumerate(zip(pos,posnext,prob)):
                if (x==xnext).all(): 
                    weight[i] += 1
                else:
                    yield mcmcsample(weight[i]+1,-l,x)
                    weight[i] = 0
            pos = posnext
