.. cosmoslik documentation master file, created by
   sphinx-quickstart on Sun Jul  1 23:05:22 2012.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

CosmoSlik Documentation
=======================

*CosmoSlik is a modular Monte Carlo sampler of cosmological likelihoods.*


Contents:

.. toctree::
    :maxdepth: 1

    Installing <install>
    The Ini File <params>
    Create your own plugins <create>
    Create your own extra-galactic foregrounds models <plugins/models.egfs>
    
.. toctree::
    :maxdepth: 3
   
    plugins
   
   


Quickstart
==========

:ref:`Install <install>` CosmoSlik, then create a file ``params.ini`` containing::

    likelihoods = wmap
    models = camb
    derivers = cosmology
    samplers = metropolis_hastings
    
    logA = 3.1 [2 4 0.1]
    ns = 1 [0 2 0.1]
    theta = 0.0104 [0 1 0.0001]
    tau = 0.09 [0 1 .01]
    ombh2 = .022 [0 1 .001]
    omch2 = .11 [0 1 .01]

    output_file = chain.dat
    samples = 100000
    
Now run::

    ./cosmoslik.py params.ini
    
and you are running a chain for the WMAP likelihood! Now add some other :ref:`plugins <plugins>` or :ref:`create <create>` your own!



Summary
=======

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


