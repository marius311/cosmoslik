from cosmoslik.chains import combine_covs
from numpy import ix_

def initialize_covariance(sampled, covs=None):
    """
    Get the covariance for the set of parameters in `sampled`
    
    
    Prepare the proposal covariance based on anything passed to
    self.proposal_cov, defaulting to the `scale` of each sampled parameter
    otherwise.
    """
    in_covs = [{k:v.scale for k,v in sampled.items() if hasattr(v,'scale')}]
    if covs is not None:
        in_covs += (covs if isinstance(covs,list) else [covs])
    names, covs = combine_covs(*in_covs)
    
    missing = [s for s in sampled if s not in names]
    if missing:
        raise ValueError("Parameters %s not in covariance and no scale given."%missing)
    
    idxs = [names.index(s) for s in sampled]
    return covs[ix_(idxs,idxs)]
