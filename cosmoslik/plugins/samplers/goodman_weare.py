from numpy import *
from collections import namedtuple
import cosmoslik.mpi as mpi
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
        
        
        #hack since emcee doesn't allow returning extra info from likelihood
        #===
        def mylnl(x):
            l, self.newp[tuple(x)] = lnl(x,p)
            return -l
        #===
        
        sampler=EnsembleSampler(nwalkers,len(x), mylnl, pool=namedtuple('pool',['map'])(mpi.mpi_map))
        p0=mpi.mpi_consistent(multivariate_normal(x,p['_cov']/len(x)*p.get('proposal_scale',1)**2,size=nwalkers))
    
        for pos,lnprob,_ in sampler.sample(p0,iterations=nsamp/nwalkers):

            #===
            if mpi.get_size()>1:
                for np in self.newp.values():
                    for k in list(np.keys()): 
                        if k[0]=='_': np.pop(k)
                newps = MPI.COMM_WORLD.gather(self.newp)
                if mpi.is_master(): 
                    for newp in newps: self.newp.update(newp) 
            self.newp = {tuple(x):self.newp.pop(tuple(x)) for x in pos} if mpi.is_master() else {}
            #===
            
            for cur_x, cur_lnl in zip(pos.copy(),lnprob.copy()):
                if mpi.is_master(): yield cur_x, 1, -cur_lnl, self.newp[tuple(cur_x)]
