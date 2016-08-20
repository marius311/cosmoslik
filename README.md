# CosmoSlik

CosmoSlik is a <b>Cosmo</b>logy <b>S</b>ampler of <b>Lik</b>elihoods. 

CosmoSlik is a simple but powerful way to quickly code up, run, and analyze a cosmology MCMC chain which can involve things like [CAMB](http://www.camb.info), [Class](http://www.class-code.net), the Planck likelihood, other cosmological likelihoods, your own customizations, etc... It has advantages in ease-of-use, flexibility, and modularity as compared to similar codes like [CosmoMC](), [MontePython](), or [CosmoSIS](). 


## Quickstart

### Basic Examples

To create a *CosmoSlik script* which runs a chain, you start with the following boilerplate code, which you should put in some file, e.g. `myscipt.py`,

```python
from cosmoslik import *

class myscript(SlikPlugin):
    def __init__(self):
        super(SlikPlugin,self).__init__()
        # initialization code here
        
    def __call__(self):
        # return negative log-likelihood here
```

This script is like the "ini" file other similar tools use to define a particular chain, but is of course much more flexible because its just Python code which can do arbitrarily complex things. 

Lets code up a simple example likelihood that samples a 2D unit Gaussian. In the initialization section we define the parameters to be sampled by the MCMC chain, and in the likelihood section we use these parameters to return a likelihood (by convention CosmoSlik wants the negative log-likelihood). Finally, we set which MCMC sampler to use. The new script looks like this:

```python
from cosmoslik import *

class myscript(SlikPlugin):
    def __init__(self):
        super(SlikPlugin,self).__init__()
        self.a = param(start=0, scale=1)
        self.b = param(start=0, scale=1)
        self.sampler = get_plugin("samplers.metropolis_hastings")(self,num_samples=1000, output_file="myscript.chain")
        
    def __call__(self):
        return (self.a**2 + self.b**2)/2
```


You can now immediately run this chain from the command line by doing, 

```bash
$ python -m cosmoslik -n 8 myscript.py
```

The `-n 8` option runs 8 chains in parallel with MPI (you will need `mpi4py` installed) and automatically updates the proposal covariance. The resulting 8 chains are all written to the single file specified in the script, in our case `myscript.chain`. It is safe to read the file while the chain is in progess, and you can do so with,

```python
>>> chain = load_chain("myscript.chain")
```

Alternatively, you could run this script from an interactive session and get a chain directly via:

```python
>>> chain = run_chain("myscript.py", nchains=8)
```

`Chain` objects of the kind returned by `load_chain` and `run_chain` have a number of useful methods for manipulating, post-processing, and plotting results. For example, since we ran 8 chains, we might want to concatenate them into a single chain, first burning-off some samples at the beginning of each chain, and then create a "triangle" plot from the result. This can be done with,

```python
>>> chain.burnin(500).join().likegrid()
```

### Script options and flexibility

The power in using Python to define scripts is that we can do lots of advanced things that can't be done in "ini" files or other mini-languages. This means a single script can in fact be a powerful multi-purpose tool that runs several different chains. For example, we can decide on the fly how many parameters to sample. Consider the following script which samples an `ndim`-dimensional Gaussian where `ndim` is a parameter we will pass in, 

```python
from cosmoslik import *

class myscript(SlikPlugin):
    def __init__(self, ndim, num_samples=1000):
        super(SlikPlugin,self).__init__(ndim=ndim)
        for i in range(ndim):
            self["param%i"%i] = param(start=0, scale=1)
        self.sampler = get_plugin("samplers.metropolis_hastings")(self, num_samples=num_samples, output_file="myscript_ndim_%i.chain"%ndim)
        
    def __call__(self):
        return sum([self["param%i"%i]**2 for i in range(self.ndim)])/2
```

A number of things are on display here:

* Giving the `__init__` function parameters, in this case the required `ndim` and the optional `num_samples` (which is then passed to the sampler).
* Defining sampled parameters dynamically from a `for` loop.
* Treating `self` as a dictionary, i.e. `self["param"]` in place of `self.param` which allows us to create the numbered parameters
* Changing the name of the `output_file` based on parameters. 

Additionally, simply by having written this script, CosmoSlik creates a nice command line wrapper you can use to specify parameters. You can see the "help" for it via:


```
$ python -m cosmoslik myscript.py -h
usage: python -m cosmoslik myscript.py [-h] [--num_samples NUM_SAMPLES] ndim

positional arguments:
  ndim

optional arguments:
  -h, --help            show this help message and exit
  --num_samples NUM_SAMPLES
                        default: 1000
```

You can then run e.g. a a 10-dimensional chain with 10000 steps via:

```bash
$ python -m cosmoslik myscript.py 10 --num_samples 10000
```

### Cosmological Examples

Having seen the basic machinery of how to write CosmoSlik scripts and run them, lets see how to run a real cosmology chain. As an example, we'll run a Planck chain, using CAMB to compute the CMB Cl's. The script file looks like this,

```python

class planck(SlikPlugin):

    def __init__(self):
        super(planck,self).__init__(self)

        self.cosmo = SlikDict(
            logA  = param(3.108,   0.03),
            ns    = param(0.962,   0.006),
            ombh2 = param(0.02221, 0.0002),
            omch2 = param(0.1203,  0.002),
            theta = param(0.0104,  0.00003),
            tau   = param(0.055,   1),
            pivot_scalar = 0.05,
            mnu = 0.06
        )
        
        self.calPlanck = param(1,0.0025,gaussian_prior=(1,0.0025))
        
        self.camb = get_plugin("models.camb")()

        self.highlTT = get_plugin('likelihoods.clik')(clik_file='plik_lite_v18_TT.clik')
        self.lowlTT = get_plugin('likelihoods.clik')(clik_file='commander_rc2_v1.1_l2_29_B.clik')
        self.priors = get_plugin('likelihoods.priors')(self)

        self.sampler = get_plugin('samplers.metropolis_hastings')(
            self,
            num_samples = 1e7,
            output_file = 'planck.chain',
       )


    def __call__(self):
        self.cosmo.As = exp(self.cosmo.logA)*1e-10
        self.highlTT.A_Planck = self.lowlTT.A_planck = self.calPlanck
        
        self.cls = self.camb(**self.cosmo)

        return lsum(
            lambda: self.priors(self),
            lambda: self.highlTT(self.cls),
            lambda: self.lowlTT(self.cls)
        )
```

Some new things here are:

* "Attaching" sampled parameters not directly to `self` but to a sub-attribute, in this case, `self.cosmo`. CosmoSlik will find all sampled parameters if they are attached to any `SlikDict`s attached to `self` (recursively, any number of `SlikDict`s deep). 
* Using the `get_plugin` function, which returns one of the many CosmoSlik plugins, e.g. the one which wraps CAMB, or which loads Planck `clik` files and uses them to evaluate the likelihood. 
* Using the `likelihoods.priors` plugin which reads through sampled parameters and imposes priors, for example the one specified by `gaussian_prior=(1,0.0025)`. 
* Using the `lsum` function. This is a convenience function which simply sums up its arguments in order, but as soon as one of them returns `inf`, it doesn't waste time evaluating the rest. 
