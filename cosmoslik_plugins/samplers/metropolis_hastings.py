from numpy import log, mean, array, sqrt, diag, genfromtxt, sum, dot, cov, inf, loadtxt, diag, nan, hstack
from numpy.random import multivariate_normal, uniform, seed
import cosmoslik.mpi as mpi
import re, time
from itertools import product
from hashlib import md5
import cPickle
from collections import defaultdict
from cosmoslik import SlikSampler, SlikFunction, param, all_kw
from cosmoslik.chains.chains import Chain, Chains
from cosmoslik.cosmoslik import sample
import multiprocessing
import struct
    
    
class mcmc_sample(sample):
    def __init__(self,weight,*args,**kwargs):
        self.weight = weight
        super(mcmc_sample,self).__init__(*args,**kwargs)


@SlikFunction
def load_chain(output_file):
    """
    Load a chain produced by metropolis_hastings2
    """
    c = cPickle.load(open(output_file,'rb'))
    if isinstance(c,Chains): return c
    else: 
        try:
            with open(output_file,'rb') as f:
                params, derived = cPickle.load(f)
                names = ['lnl','weight']+params+derived
                return Chains([Chain(zip(names,array([hstack([s.lnl,s.weight,s.x,s.extra]) for s in sample]).T)) for sample in cPickle.load(f)])
        except:
            return _load_chain_old(output_file)


@SlikFunction
def _load_chain_old(output_file):
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
        for i,k in enumerate(derived): 
            if k not in params:
                chains[source][k] += [s.extra[i] for s in samples]
        chains[source]['lnl'] += [s.lnl for s in samples]
        chains[source]['weight'] += [s.weight for s in samples]
    
    return Chains([Chain({k:array(v) for k,v in chains[i].items()}) for i in sorted(chains.keys())])

    
class metropolis_hastings(SlikSampler):
    """
    
    An adaptive metropolis hastings sampler. 
    
    To run with MPI, call your script with:

        mpiexec -n <nchains+1> python -m cosmoslik script.py
    
    (Note one process is a "master," so run one more process than you want chains)
    
    
    Args:
        params: The script to which this sampler is attached
        output_file: File where to save the chain (if running with MPI, everything still
            gets dumped into one file). The file is a Python pickle of a `Chains` object.
            Can use `cosmoslik.utils.load_chain` as a convenient way to load chains without
            having to import cPickle.
        output_extra_params: Extra parameter names besides the sampled ones which to output to file. 
        num_samples: The number of total desired samples (including rejected ones) 
        print_level: 0/1/2 to print little/medium/alot 
        proposal_cov: Path to a file which holds the proposal covariance. The format
            is a standard ascii matrix, with the first line being 
            a comment with space-separated variable names. (See e.g. savecov). This should
            be a best estimate of the posterior covariance. The actual proposal covariance is 
            multiplied by `proposal_scale**2 / N` where `N` is the number of parameters.
            (default: diagonal covariance taken from the `scale` of each parameter)
        proposal_scale: Scale the proposal matrix. (default: 2.4)
        proposal_update: Whether to update the proposal matrix. 
            Ignored if not running with MPI. The proposal is updated by taking 
            the sample covariance of the last half of each chain. (default: True) 
        proposal_update_start: If `proposal_update` is True, how many total samples (including rejected) 
            per chain to wait before starting to do updates (default: 1000). 
        mpi_comm_freq: Number of accepted samples to wait inbetween the chains
            communicating with the master process (default: 50).
        reseed: Draw a random seed based on system time and process number before starting. (default: True) 
        yield_rejected: Yield samples with 0 weigth (default: False)
        
    """
       
    
    def __init__(self, 
                 params,
                 output_file=None,
                 output_extra_params=None,
                 num_samples=100,
                 print_level=0,
                 proposal_cov=None,
                 proposal_scale=2.4,
                 proposal_update=True,
                 proposal_update_start=1000,
                 mpi_comm_freq=100,
                 temp=1,
                 reseed=True,
                 yield_rejected=False):
        """
        """
        if output_extra_params is None: output_extra_params = []
        super(metropolis_hastings,self).__init__(**all_kw(locals(),['params']))
        
        self.sampled = params.find_sampled()
        self.x0 = [params[k].start for k in self.sampled]
        self.proposal_cov = self.initialize_covariance(self.sampled)

    
    def initialize_covariance(self,sampled):
        """Load the sigma, defaulting to diagonal entries from the WIDTH of each parameter."""
        if (self.proposal_cov is None): 
            prop_names, prop = [], None
        else: 
            with open(self.proposal_cov) as f:
                prop_names = re.sub("#","",f.readline()).split()
                prop = loadtxt(f)
                
        for k,v in sampled.items():
            if not (k in prop_names or hasattr(v,'scale')): 
                raise ValueError("Parameter '%s' not in covariance and no scale given."%k)
        
        sigma = diag([getattr(v,'scale')**2. if hasattr(v,'scale') else nan for v in sampled.values()])
        common = set(sampled.keys()) & set(prop_names)
        if common:
            idxs = zip(*list(list(product([ps.index(n) for n in common],repeat=2)) for ps in [sampled.keys(),prop_names]))
            for ((i,j),(k,l)) in idxs: sigma[i,j] = prop[k,l]
            
        return sigma
  

        
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
        if mpi.is_master() and self.print_level>=1: print 'Starting MCMC chain...'
        if self.output_file is not None:
            self._output_file = open(self.output_file,"wb")
            cPickle.dump((self.sampled.keys(),self.output_extra_params),self._output_file,protocol=2)
        
        return self._mpi_mcmc(self.x0,lnl)
    
    
    def _mcmc_withprint(self,x0,lnl):
        rank = mpi.get_rank()
        for s in self._mcmc(x0, lnl):
            if self.print_level >= 2: 
                print '\033[95mChain %i:\033[0m lnl=%.2f, weight=%i, params={%s}'%\
                    (rank,s.lnl, s.weight,
                     ', '.join('%s=%.3g'%i for i in zip(self.sampled.keys(),s.x)))
            yield s
        
    def check_seed(self):
        if self.reseed: 
            seed(int(md5(str(time.time())+str(multiprocessing.current_process().pid)).hexdigest()[:8],base=16))


    def _mcmc(self,x0,lnl):

        self.check_seed()
       
        cur_lnl, cur_weight, cur_x, cur_extra = inf, 0, x0, None
        
        for _ in range(int(self.num_samples/float(max(1,mpi.get_size()-1)))):
            test_x = multivariate_normal(cur_x,self.proposal_cov/len(x0)*self.proposal_scale**2)
            test_lnl, test_extra = lnl(*test_x)
                    
            if (log(uniform()) < self.temp*(cur_lnl-test_lnl)):
                if cur_lnl!=inf: 
                    yield(mcmc_sample(cur_weight, cur_x, cur_lnl, cur_extra))
                cur_lnl, cur_weight, cur_x, cur_extra = test_lnl, 1, test_x, test_extra
            else:
                if self.yield_rejected: yield(mcmc_sample(0,test_x,test_lnl,test_extra))
                cur_weight += 1
                
    def _print_chain_stats(self,rank,samples):
        acc = sum(1 for s in samples if s.weight!=0)
        tot = sum(s.weight for s in samples)
        print '\033[1m\033[95mChain %i:\033[0m %i/%i(%.1f%%) best=%.2f last={%s}'%\
            (rank,
             acc,
             tot,
             100*float(acc)/tot,
             min(s.lnl for s in samples),
             ', '.join('%s=%.3g'%i for i in zip(self.sampled.keys(),samples[-1].x)))

    def _mpi_mcmc(self,x,lnl):  

        output_start = self._output_file.tell()
    
        (rank,size,comm) = mpi.get_mpi()
        if mpi.get_size()==1:
            sampler = self._mcmc_withprint(x, lnl)
            while True:
                samples = []
                for _ in range(self.mpi_comm_freq):
                    try:
                        s = sampler.next()
                    except Exception:
                        if self.output_file is not None:
                            self._output_file.close()
                        raise

                    yield s
                    extra=array([s.extra[k] for k in self.output_extra_params])
                    samples.append(mcmc_sample(s.weight,s.x,s.lnl,extra))
                    
                if self.output_file is not None:
                    cPickle.dump([samples],self._output_file,protocol=2)
            
        else:
        
            from mpi4py import MPI
            
            if rank==0:
                finished = [False]*(size-1)
                samples = [[] for _ in range(size-1)]
                while (not all(finished)):
                    while not comm.Iprobe(source=MPI.ANY_SOURCE, tag=0): time.sleep(.001) #Hack so OpenMPI doesn't eat 100% CPU
                    (source,new_samples)=comm.recv(source=MPI.ANY_SOURCE)
                    if (new_samples!=None): 
                        samples[source-1] += new_samples

                        t=time.time()
                        if self.proposal_update and sum(s.weight for s in samples[source-1])>self.proposal_update_start:
                            comm.send({"proposal_cov":get_new_cov(samples)},dest=source)
                        else: 
                            comm.send({},dest=source)
                        covtime = int(1e3*(time.time() - t))
                            
                        if self.print_level>=1: 
                            self._print_chain_stats(source,samples[source-1])
                            
                        if self.output_file is not None:
                            t=time.time()
                            self._output_file.seek(output_start)
                            cPickle.dump(samples,self._output_file,protocol=2)
                            self._output_file.truncate()
                            self._output_file.flush()
                            dumptime = int(1e3*(time.time() - t))
                            print '\033[93mWork for %i: propsoal=%ims dump=%ims\033[0m'%(source,covtime,dumptime) 
                    else: 
                        finished[source-1]=True
                        
            else:
                samples = []
                for s in self._mcmc_withprint(x, lnl):
                    yield s
                    extra=array([s.extra[k] for k in self.output_extra_params])
                    samples.append(mcmc_sample(s.weight,s.x,s.lnl,extra))
                    if len(samples)==self.mpi_comm_freq:
                        t = time.time()
                        comm.send((rank,samples),dest=0)
                        self.__dict__.update(comm.recv(source=0))
                        samples = []
                        print '\033[93mChain %i wasted %ims \033[0m'%(rank,int(1e3*(time.time()-t)))
                        # print sqrt(diag(self.proposal_cov))
                        
                comm.send((rank,samples),dest=0)
                comm.send((rank,None),dest=0)

       
       
def get_covariance(data,weights=None):
    if (weights==None): return cov(data.T)
    else:
        mean = sum(data.T*weights,axis=1)/sum(weights)
        zdata = data-mean
        return dot(zdata.T*weights,zdata)/(sum(weights)-1)


def get_new_cov(samples):
    data    = array(reduce(lambda a,b: a+b,[[s.x      for s in sample[len(samples)/2:]] for sample in samples]))
    weights = array(reduce(lambda a,b: a+b,[[s.weight for s in sample[len(samples)/2:]] for sample in samples]))
    return get_covariance(data, weights)
    
    
def gelman_rubin_R(samples):
    return 1
    
    
