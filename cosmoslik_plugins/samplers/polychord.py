import os.path as osp
from cosmoslik import *
from cosmoslik import mpi
from PyPolyChord import run_polychord
from PyPolyChord.settings import PolyChordSettings
from PyPolyChord.priors import UniformPrior

class polychord(SlikSampler):
    
    def __init__(
        self,
        params,
        output_file,
        nlive=None,
        num_repeats=None,
        do_clustering=None,
        read_resume=None,
        **kwargs
    ):
    
        super().__init__(**arguments(exclude=['params']))

        self.sampled = params.find_sampled()
        
        settings = self.settings = PolyChordSettings(len(self.sampled), 0)
        settings.base_dir = osp.dirname(output_file)
        settings.file_root = osp.basename(output_file)
        for k in ['nlive','num_repeats','do_clustering','read_resume']:
            v = arguments()[k]
            if v is not None:
                setattr(settings,k,v)
                 
    def sample(self, lnl):
        
        def likelihood(theta):
            return -lnl(*theta)[0], []

        def prior(hypercube):
            theta = [0.0] * len(self.sampled)
            for i,(p,x) in enumerate(zip(self.sampled.values(),hypercube)):
                if hasattr(p,"min") and hasattr(p,"max"):
                    theta[i] = UniformPrior(p.min,p.max)(x)
                else:
                    raise Exception("To use PolyChord, all parameters must have min/max specified")
            return theta

        output = run_polychord(likelihood, len(self.sampled), 0, self.settings, prior)
        output.make_paramnames_files([(k,k) for k in self.sampled.keys()])
