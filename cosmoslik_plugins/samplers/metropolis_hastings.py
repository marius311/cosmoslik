from numpy import log, mean, array, sqrt, diag, genfromtxt, sum, dot, cov, inf, loadtxt, diag, nan, hstack, empty, vstack, dtype, concatenate, ix_
from numpy.random import multivariate_normal, uniform, seed
import cosmoslik.mpi as mpi
import re, time
from itertools import product, chain
from hashlib import md5
import cPickle
from collections import defaultdict, OrderedDict
from cosmoslik import SlikSampler, SlikFunction, param, all_kw
from cosmoslik.chains.chains import Chain, Chains, combine_covs
from cosmoslik.cosmoslik import sample
import multiprocessing
import struct
    
    
class mcmc_sample(sample):
    def __init__(self,weight,*args,**kwargs):
        self.weight = weight
        super(mcmc_sample,self).__init__(*args,**kwargs)


@SlikFunction
def load_chain(filename):
    """
    Load a chain produced by the CosmoSlik metropolis_hastings sampler. 
    """
    c = cPickle.load(open(filename,'rb'))
    if isinstance(c,(Chain,Chains)): return c
    elif isinstance(c,tuple): return _load_chain_old(filename)
    else:
        with open(filename) as f:
            names = cPickle.load(f)
            dat = []
            while True:
                try:
                    dat.append(cPickle.load(f))
                except:
                    break
            ii=set(i for i,_ in dat)

            if dat[0][1].dtype.kind=='V':
                return Chains([Chain({n:concatenate([d['f%i'%k] for j,d in dat if i==j]) for k,n in enumerate(names)}) for i in ii])
            else:
                return Chains([Chain(dict(zip(names,vstack([d for j,d in dat if i==j]).T))) for i in ii])


@SlikFunction
def _load_chain_old(output_file):
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
            return _load_chain_old2(output_file)


@SlikFunction
def _load_chain_old2(output_file):
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
                 debug_output=False,
                 temp=1,
                 reseed=True,
                 yield_rejected=False):
        """
        Args:
            params: 
                The script to which this sampler is attached
            output_file: 
                File where to save the chain (if running with MPI, everything
                still gets dumped into one file). By default only sampled
                parameters get saved.  Use `cosmoslik.utils.load_chain` to load
                chains. 
            output_extra_params: 
                Extra parameters besides the sampled ones which to save to file.
                Arbitrary objects can be outputted, in which case entires should
                be tuples of (<name>,'object'), or for more efficient and faster
                file write/reads (<name>,<dtype>) where <dtype> is a valid numpy
                dtype (e.g. '(10,10)d' for a 10x10 array of doubles, etc...)
            num_samples: 
                The number of total desired samples (including rejected ones) 
            print_level: 
                0/1/2 to print little/medium/alot 
            proposal_cov: 
                One or a list of covariances which will be combined with
                K.chains.combine_cov (see documentation there for understood
                formats) to produce the full proposal covariance.  Covariance
                for any sampled parameter not provided here will be taken  from
                the `scale` attribute of that parameters. This should be a best
                estimate of the posterior covariance. The actual proposal
                covariance is  multiplied by `proposal_scale**2 / N` where `N`
                is the number of parameters. (default: diagonal covariance taken
                from the `scale` of each parameter)
            proposal_scale: 
                Scale the proposal matrix. (default: 2.4)
            proposal_update: 
                Whether to update the proposal matrix. Ignored if not running
                with MPI. The proposal is updated by taking the sample
                covariance of the last half of each chain. (default: True) 
            proposal_update_start: 
                If `proposal_update` is True, how many total samples (including
                rejected) per chain to wait before starting to do updates
                (default: 1000). 
            mpi_comm_freq:
                Number of accepted samples to wait inbetween the chains
                communicating with the master process and having their progress
                written to file (default: 50)
            reseed: 
                Draw a random seed based on system time and process number
                before starting. (default: True) 
            yield_rejected: 
                Yield samples with 0 weight (default: False)
            debug_output: 
                Print (code) debugging messages.
        """

        if output_extra_params is None: output_extra_params = []
        super(metropolis_hastings,self).__init__(**all_kw(locals(),['params']))

        self.output_extra_params = OrderedDict([k if isinstance(k,tuple) else (k,dtype('float').name) for k in output_extra_params])
        self.sampled = params.find_sampled()
        self.x0 = [params[k].start for k in self.sampled]
        self.proposal_cov = self.initialize_covariance(self.sampled)

    
    def initialize_covariance(self, sampled):
        """
        Prepare the proposal covariance based on anything passed to
        self.proposal_cov, defaulting to the `scale` of each sampled parameter
        otherwise.
        """
        in_covs = [{k:v.scale for k,v in sampled.items() if hasattr(v,'scale')}]
        if self.proposal_cov is not None:
            in_covs += (self.proposal_cov if isinstance(self.proposal_cov,list) else [self.proposal_cov])
        names, covs = combine_covs(*in_covs)
        
        missing = [s for s in sampled if s not in names]
        if missing:
            raise ValueError("Parameters %s not in covariance and no scale given."%missing)
        
        idxs = [names.index(s) for s in sampled]
        return covs[ix_(idxs,idxs)]
        
        
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
        acc = sum(1 for s in samples['f1'] if s>0)
        tot = samples['f1'].sum()
        print '\033[1m\033[95mChain %i:\033[0m %i/%i(%.1f%%) best=%.2f last={%s}'%\
            (rank,
             acc,
             tot,
             100*float(acc)/tot,
             samples['f0'].min(),
             ', '.join('%s=%.3g'%i for i in zip(self.sampled.keys(),samples[['f%i'%i for i in range(2,2+len(self.sampled))]][-1])))

    def _mpi_mcmc(self,x,lnl):  
    
        if self.output_file is not None:
            self._output_file = open(self.output_file,"wb")
            cPickle.dump(hstack([['lnl','weight'],
                                 self.sampled.keys(),
                                 self.output_extra_params.keys()]),
                         self._output_file,protocol=2)

        (rank,size,comm) = mpi.get_mpi()
        if rank==0 and size>1:
            from mpi4py import MPI
            finished = {i:False for i in range(1,size)}
            samples = {i:None for i in range(1,size)}
            while not all(finished.values()):
                while not comm.Iprobe(source=MPI.ANY_SOURCE, tag=0): time.sleep(.001) #Hack so OpenMPI doesn't eat 100% CPU
                (source,new_samples)=comm.recv(source=MPI.ANY_SOURCE)
                if (new_samples!=None): 
                    samples[source] = hstack([samples[source],new_samples]) if samples[source] is not None else new_samples

                    t=time.time()
                    if self.proposal_update and samples[source]['f1'].sum()>self.proposal_update_start:
                        comm.send({"proposal_cov":get_new_cov(samples.values(), len(self.sampled))},dest=source)
                    else: 
                        comm.send({},dest=source)
                    covtime = int(1e3*(time.time() - t))
                        
                    if self.print_level>=1: 
                        self._print_chain_stats(source,samples[source])
                        
                    if self.output_file is not None:
                        t=time.time()
                        cPickle.dump((source,new_samples),self._output_file,protocol=2)
                        self._output_file.flush()
                        dumptime = int(1e3*(time.time() - t))
                        if self.debug_output: 
                            print '\033[93mWork for %i: propsoal=%ims dump=%ims\033[0m'%(source,covtime,dumptime) 
                else: 
                    finished[source]=True
                    
        else:
            sampler = self._mcmc_withprint(x, lnl)
            dtypes = ','.join([dtype(float).name]*(2+len(self.sampled))+self.output_extra_params.values())
            samples = empty(self.mpi_comm_freq,dtypes)
            try:
                while True:
                    i=0
                    for i, s in zip(range(self.mpi_comm_freq),sampler):
                        yield s
                        samples[i] = tuple(chain(*[[s.lnl,s.weight],s.x,[s.extra[k] for k in self.output_extra_params]]))

                    if size>1:
                        t = time.time()
                        comm.send((rank,samples[:i+1]),dest=0)
                        self.__dict__.update(comm.recv(source=0))
                        if self.debug_output: 
                            print '\033[93mChain %i wasted %ims \033[0m'%(rank,int(1e3*(time.time()-t)))
                    elif self.output_file is not None:
                        cPickle.dump((0,samples[:i]),self._output_file,protocol=2)

                    if i<self.mpi_comm_freq-1: 
                        if size>1: comm.send((rank,None),dest=0)
                        break
            finally:
                if size==1 and self.output_file is not None:
                    self._output_file.close()

       
       
def get_covariance(data,weights=None):
    if (weights==None): return cov(data.T)
    else:
        mean = sum(data.T*weights,axis=1)/sum(weights)
        zdata = data-mean
        return dot(zdata.T*weights,zdata)/(sum(weights)-1)


def get_new_cov(samples, nparams):
    data = concatenate([s[['f%i'%i for i in range(2,2+nparams)]][s.shape[0]/2:] for s in samples if s is not None])
    data = data.view(float).reshape(data.shape+(-1,))
    weights = hstack([s['f1'][s.shape[0]/2:] for s in samples if s is not None])
    return get_covariance(data, weights)
    
def gelman_rubin_R(samples):
    return 1
    
    
