from numpy import log, mean, array, sqrt, diag, genfromtxt, sum, dot, cov, inf
from random import random, choice
from numpy.random import multivariate_normal
import cosmoslik.mpi as mpi
import re, time
from itertools import product
from collections import namedtuple
from cosmoslik.plugins import Sampler
from cosmoslik.params import SectionDict
    
sampletuple = namedtuple('sampletuple',['x','weight','lnl','extra'])
    

class bpmc(Sampler):
    """
    
    =============================
    Branch Prediction Monte-Carlo
    =============================
    
    Usage
    =====
    
    To use this plugin set ``samplers = bpmc`` 
    
    
    Description
    ===========
    
    If you run multiple MPI processes, they work together on a single Metropolis-Hastings chain
    using branch prediction. As long as the typical weight per step is more than the number 
    of processes, you'll get linear speed-up, even during burn-in. 
    
    
    Parameters
    ==========
    
    proposal_matrix
    ---------------
        Path to a file which holds the proposal covariance. The format
        is a standard ascii matrix, with the first line being 
        a comment with space-separated variable names. (See e.g. savecov). This should
        be a best estimate of the posterior covariance. The actual proposal covariance is 
        multiplied by ``proposal_scale**2 / N`` where ``N`` is the number of parameters.
        (default: diagonal covariance taken from the ``width`` of each parameter)
        
    proposal_scale
    --------------
        Scale the proposal matrix. (default: ``2.4``)
        
    proposal_update
    ---------------
        Whether to update the proposal based on the sample covariance of previous steps. 
        Ignored if not running with MPI. The proposal is updated by taking 
        the sample covariance of the last half of each chain. (default: ``True``) 
        
    proposal_update_start
    ---------------------
        If ``proposal_update`` is ``True``, how many non-unique samples per chain to
        wait before starting to do updates. (default: ``1000``) 
        
    
        
    
    
    """

    def init(self,p):
        if mpi.get_size()>1 and mpi.get_rank()>0: p.pop('output_file')
        
    def sample(self,x,lnl,p):
        """
        Returns a generator which yields samples from the likelihood 
        using the Metropolis-Hastings algorithm.
        
        The samples returned are tuples of (x,weight,lnl,extra)
            lnl - likelihood
            weight - the statistical weight (could be 0 for rejected steps)
            x - the vector of parameter values
            extra - the extra information returned by lnl
        """
        
        return _mcmc(x,lnl,p)
    
    
def _mcmc(x,lnl,p):
   
    cur_lnl, cur_weight, cur_x, cur_extra = inf, 0, x, None

    samples, weights = [], [] 
    
    def update(v): 
        if v!=None: p.update(v)
        
    def makepicklable(extra):
        return {} if not extra else SectionDict({k:v for k,v in extra.items() if k[0]!='_'})

    def step((cur_lnl, test_x)):
        test_lnl, test_extra = lnl(test_x, p)
        return (log(random()) < cur_lnl-test_lnl), sampletuple(x=test_x, weight=None, lnl=test_lnl, extra=makepicklable(test_extra))


    for i in range(p.get("samples",100)):
        
        cur_lnl, cur_x, cur_extra = mpi.mpi_consistent((cur_lnl, cur_x, makepicklable(cur_extra)))
        
        test_x = [multivariate_normal(cur_x,p['_cov']/len(x)*p.get('proposal_scale',2.4)**2) for _ in range(mpi.get_size())]
        branches = mpi.mpi_map(step, zip([cur_lnl]*mpi.get_size(),test_x), distribute=False)

        if mpi.is_master():
            if i%500==0 and i>1000: 
                p['_cov'] = get_new_cov(samples, weights)
            for accepted, sample in branches:
                if not accepted:
                    update((yield(sampletuple(x=sample.x, weight=0, lnl=sample.lnl, extra = sample.extra))))

            cur_weight += sum(1 for accepted,_ in branches if not accepted)

            if any(accepted for accepted,_ in branches):
                if cur_lnl!=inf: 
                    update((yield(sampletuple(x=cur_x, weight=cur_weight, lnl=cur_lnl, extra=cur_extra))))
                    samples.append(cur_x)
                    weights.append(cur_weight)
                cur = choice([sample for accepted,sample in branches if accepted])
                cur_lnl, cur_x, cur_extra = cur.lnl, cur.x, cur.extra 
                cur_weight = 0
        
            
                  
def get_covariance(data,weights=None):
    if (weights==None): return cov(data.T)
    else:
        mean = sum(data.T*weights,axis=1)/sum(weights)
        zdata = data-mean
        return dot(zdata.T*weights,zdata)/(sum(weights)-1)


def get_new_cov(samples,weights=None):
    """
    shape(samples) = (nchains,nsamples,nparams)
    shape(weights) = (nchains,nsamples)
    """
    n = len(samples)/2
    return get_covariance(array(samples[n:]), weights[n:])
    
    
def gelman_rubin_R(samples):
    return 1
    
