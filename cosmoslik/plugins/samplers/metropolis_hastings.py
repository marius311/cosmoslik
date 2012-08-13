from numpy import log, mean, array, sqrt, diag, genfromtxt, sum, dot, cov, inf
from random import random
from numpy.random import multivariate_normal
import cosmoslik.mpi as mpi, re, time
from itertools import product
from collections import namedtuple
from cosmoslik.plugins import Sampler
    
sampletuple = namedtuple('sampletuple',['x','weight','lnl','extra'])
    

class metropolis_hastings(Sampler):
    """
    
    ===================
    Metropolis-Hastings
    ===================
    
    
    Usage
    =====
    
    To use this module set ``samplers = metropolis_hastings`` 
    
    
    Running with MPI
    ================
    
    This sampler can be run with MPI. For example to run 8 parallel chains use::
    
        python -m cosmoslik -n 8 params.ini
    
    or::
    
        mpiexec -n 9 python -m cosmoslik params.ini
    
    (When running ``mpiexec`` by hand, one process is a "master," so run one extra process.)
    
    
    
    Parameters
    ==========
    
    include sampler ones...
    
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
        if mpi.get_rank()>0: 
            if 'output_file' in p: p['output_file']+=('_%i'%mpi.get_rank())
            if 'samples' in p: p['samples']/=(mpi.get_size()-1)
            p['quiet']=True 
        elif mpi.get_size()>1: 
            p.pop('output_file')
        
        
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
        
        if mpi.get_size()==1: return _mcmc(x,lnl,p)
        else: return _mpi_mcmc(x,lnl,p)
    
    
    
def _mcmc(x,lnl,p):
   
    cur_lnl, cur_weight, cur_x = inf, 0, x

    def update(v): 
        if v!=None: p.update(v)
        
    for _ in range(p.get("samples",100)):
        test_x = multivariate_normal(cur_x,p['_cov']/len(x)*p.get('proposal_scale',2.4)**2)
        test_lnl, test_extra = lnl(test_x, p)
                
        if (log(random()) < cur_lnl-test_lnl):
            if cur_lnl!=inf: update((yield(sampletuple(cur_x, cur_weight, cur_lnl, cur_extra))))
            cur_lnl, cur_weight, cur_x, cur_extra = test_lnl, 1, test_x, test_extra
        else:
            if p.get('metropolis_hastings',{}).get('weighted_samples',True): cur_weight+=1
            else: update((yield(sampletuple(cur_x, cur_weight, cur_lnl, cur_extra))))
            
            if p.get('metropolis_hastings',{}).get('rejected_samples',True): 
                update((yield(sampletuple(test_x, 0, test_lnl, test_extra))))
            
            
            
def _mpi_mcmc(x,lnl,p):  

    (rank,size,comm) = mpi.get_mpi()
    from mpi4py import MPI #FIX THIS
    
    if rank==0:
        finished = [False]*(size-1)
        samples, weights, lnls = [[[] for _ in range(size-1)] for __ in range(3)] 
        while (not all(finished)):
            while not comm.Iprobe(source=MPI.ANY_SOURCE, tag=0): time.sleep(.1) #Hack so OpenMPI doesn't eat 100% CPU
            (source,new_samples)=comm.recv(source=MPI.ANY_SOURCE)
            if (new_samples!=None): 
                lnls[source-1]+=[s.lnl for s in new_samples]
                samples[source-1]+=[s.x for s in new_samples]
                weights[source-1]+=[s.weight for s in new_samples]
                
                if (p.get("proposal_update",True) 
                    and sum(weights[source-1])>p.get('proposal_update_start',1000)):
                    comm.send({"_cov":get_new_cov(samples,weights)},dest=source)
                else: comm.send({},dest=source)
                comm.send(None,dest=source)

            else: 
                finished[source-1]=True
                
    else:
        samples = []
        sampler = _mcmc(x, lnl, p)
        s=sampler.next()
        while True:
            try:
                yield s
                samples.append(sampletuple(s.x,s.weight,s.lnl,None))
                if len(samples)==50:
                    comm.send((rank,samples),dest=0)
                    s = sampler.send(comm.recv(source=0))
                    samples = []
                else: s = sampler.next()
            except StopIteration: break
        comm.send((rank,None),dest=0)

       
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
    data = array(reduce(lambda a,b: a+b,[s[len(s)/2:] for s in samples]))
    weights = array(reduce(lambda a,b: a+b,[w[len(w)/2:] for w in weights]))
    return get_covariance(data, weights)
    
    
def gelman_rubin_R(samples):
    return 1
    
