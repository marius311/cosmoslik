from pywmap import pywmap
from numpy import zeros
from cosmoslik.plugins import Likelihood
import os


class wmap(Likelihood):
    """
    ===============
    WMAP Likelihood
    ===============
    
    - Written by WMAP team (see `<http://lambda.gsfc.nasa.gov/>`_)
    - CosmoSlik module by Marius Millea
    - Updated July 1, 2012
    
    Description
    ===========
        
    This module wraps the official WMAP likelihood code. 
    Some minor modifications were made to allow:
    
    - Choosing the WMAP data directory at runtime
    - Choosing the lmin/lmax at runtime
    
    
    Install Notes
    =============
    
    This build this module run::
    
        ./cosmoslik.py --build likelihoods.wmap
        
    The Makefile for this module reads the following flags from ``Makefile.inc``:
    
    - ``$(CFITSIO)``
    - ``$(LAPACK)``
    - ``$(F2PYFLAGS)``
    
    
    Models
    ======
    
    The WMAP module requires a `Model` which provides the following:
    
    - ``cl_TT``
    - ``cl_TE``
    - ``cl_EE``
    - ``cl_BB``
        
    Extra-galactic foregrounds are ignored. 
    
    
    Parameters
    ==========
    
    This module reads the following parameters from the ini file:
    
    [wmap].data_dir
    ---------------
        The path to the wmap/data directory.
    
    [wmap].use
    ----------
        A subset of ``['TT','TE','EE','BB']`` corresponding to 
        which likelihood terms to use.
        
    [wmap].TT.lrange
    ----------------
        The TT range in ell to use in the likelihood
    
    [wmap].TE.lrange
    ----------------
        The TE range in ell to use in the likelihood
    
    """

    def init(self,p):
        self.use = p.get(('wmap','use'),['TT','TE','EE','BB'])
        ttmin, ttmax = p.get(('wmap','TT.lrange'),(2,1200))
        temin, temax = p.get(('wmap','TE.lrange'),(2,800))
        datadir = p.get(('wmap','data_dir'),None)
        if not datadir: raise Exception('Please specify the WMAP data directory in the parameter file with:\n[wmap]{\n  data_dir=/path/to/wmap/data\n}')
        elif not os.path.exists(datadir): raise Exception("The WMAP data directory you specified does not exist: '%s'"%datadir)
        pywmap.wmapinit(ttmin,ttmax,temin,temax,os.path.normpath(datadir)+'/')
    
    def get_required_models(self, p):
        return ['cl_%s'%x for x in self.use]
    
    def lnl(self, p, model):
        cltt, clte, clee, clbb = [zeros(1202) for _ in range(4)]
        
        for cl,x in zip([cltt,clte,clee,clbb],['TT','TE','EE','BB']):
            if x in self.use:
                m = model['cl_%s'%x]
                s = slice(0,min(len(m),len(cl)))
                cl[s] = m[s]
    
        liketerms = pywmap.wmaplnlike(cltt=cltt[2:],clte=clte[2:],clee=clee[2:],clbb=clbb[2:])
        return sum(liketerms)
    
    
