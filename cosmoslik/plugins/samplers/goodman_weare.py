from numpy import *
from collections import namedtuple
import cosmoslik.mpi as mpi
import cosmoslik.chains as chains
from emcee.ensemble import EnsembleSampler
from numpy.random import multivariate_normal
from cosmoslik.plugins import Sampler
from mpi4py import MPI


class goodman_weare(Sampler):

    newp = {}
    
    def sample(self,x,lnl,p):
        nwalkers = p.get(('goodman_weare','walkers'),max(1,2*len(x)))
        if nwalkers < 2*len(x): raise Exception("You need at least %i Goodman-Weare walkers (twice the number of sampled parameters)"%(2*len(x)))
        nsamp = p.get('samples',10000)
        
        def mylnl(x):
            l, extra = lnl(x,p)
            return -l, extra
        
        sampler=EnsembleSampler(nwalkers,len(x), mylnl, pool=namedtuple('pool',['map'])(mpi.mpi_map))
        
        p0=mpi.mpi_consistent(multivariate_normal(x,p['_cov']/len(x)*p.get('proposal_scale',1)**2,size=nwalkers))
        
        if ('goodman_weare','starting_points') in p:
            chain = chains.load_chain(p['goodman_weare','starting_points']).join()
            for i,k in enumerate(p.get_all_sampled()):
                name = '.'.join(k)
                if name in chain: p0[:,i] = chain[name][-nwalkers:]
    
        for pos,lnprob,_,extra in sampler.sample(p0,iterations=nsamp/nwalkers):

            for cur_x, cur_lnl, cur_extra in zip(pos,lnprob,extra):
#                import ipdb; ipdb.set_trace()
                if mpi.is_master(): yield cur_x, 1, -cur_lnl, cur_extra
