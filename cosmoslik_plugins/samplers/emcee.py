from cosmoslik import *
from cosmoslik import mpi
from cosmoslik.chains.chains import combine_covs
from emcee import EnsembleSampler
from collections import OrderedDict, namedtuple
from numpy import ix_, zeros, ceil, hstack, dtype, empty, ones
from numpy.random import multivariate_normal
import pickle
from itertools import chain

mcmcsample = namedtuple("mcmcsample",["weight","lnl","x"])

class emcee(SlikSampler):
    
    def __init__(
        self,
        params,
        num_samples,
        nwalkers=100,
        proposal_cov=None,
        output_extra_params=None,
        output_file=None,
        output_freq=10
    ):
    
        if output_extra_params is None: output_extra_params = []
        super().__init__(**arguments(exclude=['params']))

        self.output_extra_params = OrderedDict([k if isinstance(k,tuple) else (k,dtype('float').name) for k in output_extra_params])
        self.sampled = params.find_sampled()
        self.x0 = [params[k].start for k in self.sampled]
        self.proposal_cov = self.initialize_covariance(self.sampled)
        self.pool = mpi.get_pool() if mpi.get_size()>1 else None
    
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
        
        # set up output file
        if self.output_file and mpi.is_master():
            _output_file = open(self.output_file,"wb")
            protocol = pickle.HIGHEST_PROTOCOL
            pickle.dump(
                obj=list(chain(['lnl','weight'], self.sampled.keys(), self.output_extra_params.keys())),
                file=_output_file,
                protocol=protocol
            )
            dtypes = ','.join([dtype(float).name]*(2+len(self.sampled))+list(self.output_extra_params.values()))
            samples = {i:empty(self.output_freq,dtypes) for i in range(self.nwalkers)}
        
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
                        
        # step each walker once, yield them all one-by-one, repeat
        weight = ones(self.nwalkers)
        isample = zeros(self.nwalkers,dtype=int)
        sampler = EnsembleSampler(self.nwalkers, len(self.sampled), lnprob, pool=self.pool)

        nsteps = int(ceil(self.num_samples/self.nwalkers))
        for i in range(nsteps):
            posnext, prob, state, blobs = sampler.run_mcmc(pos,1)
            
            for iwalker, (x, xnext, l, params) in enumerate(zip(pos,posnext,prob,blobs)):
                if (x==xnext).all() and i!=nsteps-1: 
                    weight[iwalker] += 1
                else:
                    yield mcmcsample(weight[iwalker],-l,x)
                    
                    # write to file once every `self.output_freq` accepted steps (per walker)
                    if self.output_file and mpi.is_master():
                        row = tuple(chain([-l,weight[iwalker]],x,[params[k] for k in self.output_extra_params]))
                        samples[iwalker][isample[iwalker]] = row
                        isample[iwalker] += 1
                        if isample[iwalker]>=self.output_freq or i==nsteps-1:
                            pickle.dump((iwalker,samples[iwalker][:isample[iwalker]]),_output_file,protocol)
                            isample[iwalker] = 0
                    
                    weight[iwalker] = 1
                    
            pos = posnext
