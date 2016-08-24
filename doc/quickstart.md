# Quickstart

## Basic Examples

To create a CosmoSlik *script* which runs a chain, you start with the following boilerplate code, which you should put in some file, e.g. `myscipt.py`,

```python
from cosmoslik import *

class myscript(SlikPlugin):
    def __init__(self):
        super().__init__()
        # your initialization code will go here
        
    def __call__(self):
        # you'll return negative log-likelihood here
```

This script is like the "ini" file other similar codes use to define a particular run, but is much more flexible because its just Python code which can do arbitrarily complex things, as we'll see below. 

Lets code up a simple example that samples a 2D unit Gaussian likelihood. In the initialization section we define the parameters to be sampled by the MCMC chain, and in the likelihood section we use these parameters to return a likelihood (by convention, CosmoSlik wants the negative log-likelihood). Finally, we set which MCMC sampler to use. The new script looks like this:

```python
from cosmoslik import *

class myscript(SlikPlugin):
    def __init__(self):
        super().__init__()
        
        # define sampled parameters
        self.a = param(start=0, scale=1)
        self.b = param(start=0, scale=1)
        
        # set the sampler
        self.sampler = samplers.metropolis_hastings(
            self,
            num_samples=1e5, 
            output_file="myscript.chain",
        )

    def __call__(self):
        # compute the likelihood
        return (self.a**2 + self.b**2)/2
```


You can now immediately run this chain from the command line by doing, 

```bash
$ cosmoslik -n 8 myscript.py
```

The `-n 8` option runs 8 chains in parallel with MPI (you will need `mpi4py` installed) and automatically updates the proposal covariance. The resulting 8 chains are all written to the single file specified in the script, in our case `myscript.chain`. You can read the chain file (including while the chain is running to see its progress) via

```python
>>> chain = load_chain("myscript.chain")
```

Alternatively, you could run this script from an interactive session and get a chain directly via

```python
>>> chain = run_chain("myscript.py", nchains=8)
```

Or you could have skipped the script file entirely and just run

```python
>>> chain = run_chain(myscript, nchains=8)
```

assuming you defined `myscript` in your interactive session. 


`load_chain` and `run_chain` return `Chain` or `Chains` objects which have a number of useful methods for manipulating, post-processing, and plotting results. For example, since we ran 8 chains, we might want to concatenate them into a single chain, first burning-off some samples at the beginning of each chain, and then create a "triangle" plot from the result. This would be done with,

```python
>>> chain.burnin(500).join().likegrid()
```

## Script options and flexibility

The power in using Python to define scripts is that we can do lots of advanced things that can't be done in "ini" files or other mini-languages. This means a single script can in fact be a powerful multi-purpose tool that runs several different chains. For example, we can decide on the fly how many parameters to sample. Consider the following script which samples an `ndim`-dimensional Gaussian where `ndim` is a parameter we will pass in, 

```python
from cosmoslik import *

class myscript(SlikPlugin):
    def __init__(self, ndim, num_samples=1000):
        super().__init__()
        # save the ndim parameter so we can use it below in __call__ too
        ndim = self.ndim = int(ndim)
        
        # dynamically create the sampled parameters we need
        for i in range(ndim):
            self["param%i"%i] = param(start=0, scale=1)
            
        self.sampler = samplers.metropolis_hastings(
            self, 
            num_samples=num_samples, 
            output_file="myscript_ndim_%i.chain"%ndim
        )
        
    def __call__(self):
        return sum([self["param%i"%i]**2 for i in range(self.ndim)])/2
```

Let's break down some new things we did here:

* We gave the `__init__` function some arguments, `ndim` and `num_samples`
* We defined the sampled parameters for this chain dynamically with a `for` loop.
* We used `self` as a dictionary, i.e. `self["param"]` in place of `self.param`, which allows us to create the numbered parameters

Additionally, simply by having written this script, CosmoSlik creates a nice wrapper you can use to call this script from the command line and specify parameters. You can see the "help" for it via:


```
$ cosmoslik myscript.py -h
usage: cosmoslik myscript.py [-h] [--num_samples NUM_SAMPLES] ndim

positional arguments:
  ndim

optional arguments:
  -h, --help            show this help message and exit
  --num_samples NUM_SAMPLES
                        default: 1000
```

You can then run e.g. a a 10-dimensional chain with 10000 steps via:

```
$ cosmoslik myscript.py 10 --num_samples 10000
```

## Cosmological Examples

Having seen the basic machinery of how to write CosmoSlik scripts and run them, lets see how to run a real cosmology chain. As an example, we'll run a Planck chain, using CAMB to compute the CMB Cl's. The script file looks like this,

```python
class planck(SlikPlugin):

    def __init__(self):
        super().__init__()

        # define sampled cosmological parameters
        param = param_shortcut('start','scale')
        self.cosmo = SlikDict(
            logA  = param(3.108,   0.03),
            ns    = param(0.962,   0.006),
            ombh2 = param(0.02221, 0.0002),
            omch2 = param(0.1203,  0.002),
            theta = param(0.0104,  0.00003),
            tau   = param(0.055,   0.01),
        )
        
        # sample the Planck calibration as well
        self.calPlanck = param(1,0.0025,gaussian_prior=(1,0.0025))
        
        # load CAMB to compute Cls
        self.camb = models.camb()

        # load the Planck likelihood files
        self.highlTT = likelihoods.clik(clik_file='plik_lite_v18_TT.clik')
        self.lowlTT = likelihoods.clik(clik_file='commander_rc2_v1.1_l2_29_B.clik')

        self.sampler = samplers.metropolis_hastings(
            self,
            num_samples = 1e7,
            output_file = 'planck.chain',
        )


    def __call__(self):
        # we sample logA but CAMB needs As
        self.cosmo.As = exp(self.cosmo.logA)*1e-10
        
        # compute Cls
        self.cls = self.camb(**self.cosmo)
        
        # the two Planck likelihoods read the calibration from `A_Planck` and
        # `A_planck`, so set those based on our sampled parameter
        self.highlTT.A_Planck = self.lowlTT.A_planck = self.calPlanck

        # compute likelihood
        return lsum(
            lambda: self.highlTT(self.cls),
            lambda: self.lowlTT(self.cls)
        )
```

Some new things here are:

* "Attaching" sampled parameters not directly to `self` but to a sub-attribute, in this case, `self.cosmo`. CosmoSlik will find all sampled parameters if they are attached to any `SlikDict`s attached to `self` (recursively, any number of `SlikDict`s deep). You can use this to organize parameters into convenient subgroups. 
* Using the `lsum` function. This is a convenience function which simply sums up its arguments in order, but if one of them returns `inf`, it doesn't waste time evaluating the rest. 
