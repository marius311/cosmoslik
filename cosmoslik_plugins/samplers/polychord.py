import os, os.path as osp
from cosmoslik import *
from cosmoslik import mpi
from PyPolyChord import run_polychord
from PyPolyChord.settings import PolyChordSettings
from PyPolyChord.priors import UniformPrior
from numpy import nan

class polychord(SlikSampler):
    
    def __init__(
        self,
        params,
        output_file,
        output_extra_params=[],
        nlive=None,
        num_repeats=None,
        do_clustering=None,
        read_resume=None,
        **kwargs
    ):
    
        super().__init__(**arguments(exclude=['params']))

        self.sampled = params.find_sampled()
        
        settings = self.settings = PolyChordSettings(len(self.sampled), len(self.output_extra_params))
        output_file = './'+osp.normpath(output_file)
        settings.base_dir = osp.dirname(output_file)
        settings.file_root = osp.basename(output_file)
        if mpi.is_master() and not osp.exists(self.settings.base_dir): os.makedirs(self.settings.base_dir)
        for k in ['nlive','num_repeats','do_clustering','read_resume']:
            v = arguments()[k]
            if v is not None:
                setattr(settings,k,v)
                 
    def sample(self, lnl):
        
        def likelihood(theta):
            l,p = lnl(*theta)
            return -l, [p.get(k,nan) for k in self.output_extra_params]

        def prior(hypercube):
            theta = [0.0] * len(self.sampled)
            for i,(p,x) in enumerate(zip(self.sampled.values(),hypercube)):
                if hasattr(p,"min") and hasattr(p,"max"):
                    theta[i] = UniformPrior(p.min,p.max)(x)
                elif hasattr(p,"range"):
                    theta[i] = UniformPrior(*p.range)(x)
                else:
                    raise Exception("To use PolyChord, all parameters must have min/max specified")
            return theta

        if mpi.is_master():
            with open(self.output_file+'.paramnames','w') as f:
                for k in list(self.sampled.keys())+self.output_extra_params:
                    f.write(k+' '+k+'\n')
                                       
        output = run_polychord(likelihood, len(self.sampled), len(self.output_extra_params), self.settings, prior)
        
        return []
