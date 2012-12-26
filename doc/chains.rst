.. _chains:

======
Chains
======

CosmoSlik contains utilities to load, plot, and post-process MCMC chains. Here's 
an example of typical usage, ::

    import cosmoslik as K
    
    #load chains from a CosmoSlik run
    chains = K.chains.load_chain('cosmoslik.ini')
    
    #view trace plot
    chains.plot('param')
    
    #join chains into one
    chain = chains.burnin(2000).join()
    
    #view parameter grid
    chain.likegrid()       
    
    #calculate mean/std-dev
    zip(chain.params(),chain.mean(),chain.std())
    
Below is more detail on all available functions.

Loading Chains
--------------

To load a chain, use,

.. autofunction:: cosmoslik.chains.load_chain


Chains Object
-------------

The result of :func:`~cosmoslik.chains.load_chain` for multiple chains 
is an object of type :class:`~cosmoslik.chains.Chains`

.. autoclass:: cosmoslik.chains.Chains
    :members:
    :undoc-members:


Chain Object
------------

The result of :func:`~cosmoslik.chains.load_chain` for a single chain, 
or the result of calling :meth:`~cosmoslik.chains.Chains.join` 
is an object of type :class:`~cosmoslik.chains.Chain`

.. autoclass:: cosmoslik.chains.Chain
    :members:
    :undoc-members:

Plotting
--------

Generally these functions aren't used directly, but are called from 
the methods in :class:`~cosmoslik.chains.Chain`. 
One exception is :func:`~cosmoslik.chains.likegrid` which can be 
used to simultaneously plot multiple chains on a single parameter grid.

.. autofunction:: cosmoslik.chains.like1d
.. autofunction:: cosmoslik.chains.like2d
.. autofunction:: cosmoslik.chains.likegrid

