from numpy import log, mean, array, sqrt, diag, genfromtxt, sum, dot, cov, inf, loadtxt, diag, nan
from random import random
from numpy.random import multivariate_normal
import cosmoslik.mpi as mpi, re, time
from itertools import product

import cPickle
from collections import defaultdict
from cosmoslik import SlikSampler, SlikFunction, param
from cosmoslik.chains.chains import Chain, Chains
from cosmoslik.cosmoslik import sample
    
    
    
class mcmc_sample(sample):
    def __init__(self,weight,*args,**kwargs):
        self.weight = weight
        super(mcmc_sample,self).__init__(*args,**kwargs)


@SlikFunction
def load_chain(output_file):
    """
    Load a chain produced by metropolis_hastings2
    """
    dat=[]
    with open(output_file) as f:
        while True:
            try: dat.append(cPickle.load(f))
            except: break
                
    params, derived = dat[0]
    chains = defaultdict(lambda: {k:[] for k in params+derived+['weight','lnl']})
    
    for source,samples in dat[1:]:
        for i,k in enumerate(params): chains[source][k] += [s.x[i] for s in samples]
        for i,k in enumerate(derived): chains[source][k] += [s.extra[i] for s in samples]
        chains[source]['lnl'] += [s.lnl for s in samples]
        chains[source]['weight'] += [s.weight for s in samples]
    
    return Chains([Chain({k:array(v) for k,v in chains[i].items()}) for i in sorted(chains.keys())])

    
class metropolis_hastings(SlikSampler):
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
       
    
    def __init__(self, 
                 output_file=None,
                 output_extra_params=None,
                 num_samples=100,
                 proposal_cov=None,
                 proposal_scale=2.4,
                 proposal_update=True,
                 proposal_update_start=1000,
                 mpi_comm_freq=50):
        """
        """
        if output_extra_params is None: output_extra_params = []
        ls = locals()
        ls.pop('self')
        super(metropolis_hastings,self).__init__(**ls)
        
    
    def initialize_covariance(self,slik):
        """Load the sigma, defaulting to diagonal entries from the WIDTH of each parameter."""
        if (self.proposal_cov is None): 
            prop_names, prop = [], None
        else: 
            with open(self.proposal_cov) as f:
                prop_names = re.sub("#","",f.readline()).split()
                prop = loadtxt(f)
                
        sampled = slik.get_sampled()
        
        for k,v in sampled.items():
            if not (k in prop_names or hasattr(v,'scale')): 
                raise ValueError("Parameter '%s' not in covariance and no scale given."%k)
        
        sigma = diag([getattr(v,'scale')**2 if hasattr(v,'scale') else nan for v in sampled.values()])
        common = set(sampled.keys()) & set(prop_names)
        if common:
            idxs = zip(*list(list(product([ps.index(n) for n in common],repeat=2)) for ps in [sampled.keys(),prop_names]))
            for ((i,j),(k,l)) in idxs: sigma[i,j] = prop[k,l]
            
        return sigma

    
    def init(self,slik):
        self.chain_metadata = (slik.get_sampled().keys(),self.output_extra_params)
        self.x0 = [slik.params[k].start for k in slik.get_sampled()]
        self.proposal_cov = self.initialize_covariance(slik)
            

        
        
    @SlikFunction
    def run(self,slikself):
        samps = list(slikself.sample())
        chain = Chain()
        for i,k in enumerate(slikself.get_sampled()):
            chain[k] = array([s.x[i] for s in samps])
        chain['weight'] = array([s.weight for s in samps])
        chain['lnl'] = array([s.lnl for s in samps])
        return chain
    

    def print_sample(self):
        pass
        #        if nsamp%getattr(self,'update_frequency',1)==0:
#            print "%saccepted=%s/%i(%.1f%%) best=%.2f last={%s}" % \
#                ('' if mpi.get_rank()==0 else 'Chain %i: '%mpi.get_rank(),
#                 len(samples.weight),
#                 nsamp,
#                 100*float(len(samples.weight))/nsamp,
#                 min(samples.lnl+[inf]),
#                 ', '.join([('like:%.2f'%l1)]+['%s:%s'%('.'.join(name),'%.4g'%self[name] if isinstance(self[name],float) else self[name])
#                                               for name in outputted])
#                 )

        
    def sample(self,lnl):
        """
        Returns a generator which yields samples from the likelihood 
        using the Metropolis-Hastings algorithm.
        
        The samples returned are tuples of (x,weight,lnl,extra)
            lnl - likelihood
            weight - the statistical weight (could be 0 for rejected steps)
            x - the vector of parameter values
            extra - the extra information returned by lnl
        """
        
        self._output_file = open(self.output_file,"w")
        cPickle.dump(self.chain_metadata,self._output_file)
        self._output_file.flush()
        
        if mpi.get_size()==1: return self._mcmc(self.x0,lnl)
        else: return self._mpi_mcmc(self.x0,lnl)
    
    
    
    def _mcmc(self,x0,lnl):
       
        cur_lnl, cur_weight, cur_x, cur_extra = inf, 0, x0, None
        
        for _ in range(int(self.num_samples/float(max(1,mpi.get_size()-1)))):
            test_x = multivariate_normal(cur_x,self.proposal_cov/len(x0)*self.proposal_scale**2)
            test_lnl, test_extra = lnl(*test_x)
                    
            if (log(random()) < cur_lnl-test_lnl):
                if cur_lnl!=inf: 
                    yield(mcmc_sample(cur_weight, cur_x, cur_lnl, cur_extra))
                    #TODO: print
                cur_lnl, cur_weight, cur_x, cur_extra = test_lnl, 1, test_x, test_extra
            else:
                cur_weight += 1
                
            
            
    def _mpi_mcmc(self,x,lnl):  
    
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
                    
                    if self.proposal_update and sum(weights[source-1])>self.proposal_update_start:
                        comm.send({"proposal_cov":get_new_cov(samples,weights)},dest=source)
                    else: 
                        comm.send({},dest=source)
                    
                    cPickle.dump((source,[s for s in new_samples if s.weight>0]),self._output_file,protocol=2)
                    self._output_file.flush()
                else: 
                    finished[source-1]=True
                    
        else:
            samples = []
            for s in self._mcmc(x, lnl):
                yield s
                extra=array([s.extra[k] for k in self.output_extra_params])
                samples.append(mcmc_sample(s.weight,s.x,s.lnl,extra))
                if len(samples)==self.mpi_comm_freq:
                    comm.send((rank,samples),dest=0)
                    self.__dict__.update(comm.recv(source=0))
                    samples = []
                    
            comm.send((rank,samples),dest=0)
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
    
    
