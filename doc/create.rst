.. _create:

===============================
How to Create CosmoSlik Plugins
===============================

CosmoSlik plugins come in four types, :class:`~cosmoslik.plugins.Likelihood`, 
:class:`~cosmoslik.plugins.Model`, :class:`~cosmoslik.plugins.Deriver`,
and :class:`~cosmoslik.plugins.Sampler`. When CosmoSlik sees a line in the
ini file such as, ::

    likelihoods = my_lnl
    
it executes the following Python command, ::

    from cosmoslik.plugins.likelihoods.my_lnl import my_lnl
    
In this case, the imported attribute ``my_lnl`` must be of type 
:class:`~cosmoslik.plugins.Likelihood`. 

So, to create your own plugin, all you need to do is,

* Add the relevant line to the ini file
* Ensure that the corresponding import statement succeeds
* Code functionality into your plugin

Ensuring that the import statement succeeds is as easy as creating a few files. 
Continuing the above example, you would create the following folder structure, ::

    |cosmoslik
    |  cosmoslik
    |    plugins
    |      likelihoods
    |        my_lnl
    |          __init__.py
    |          my_lnl.py

The file ``my_lnl.py`` should contain, ::

    from cosmoslik.plugins import Likelihood
    
    class my_lnl(Likelihood):
    
        def lnl(self, p, model):
            #my likelihood code here
            return inf
    
            
And the file ``__init__py`` should contain, ::

    from my_lnl import my_lnl


Types of Plugins
----------------

The base class for all CosmoSlik plugins is,

.. autoclass:: cosmoslik.plugins.CosmoSlikPlugin
    :members:

There are four types of plugins,

.. autoclass:: cosmoslik.plugins.Likelihood
    :members:
    
.. autoclass:: cosmoslik.plugins.Model  
    :members:
    
.. autoclass:: cosmoslik.plugins.Deriver
    :members:
    
.. autoclass:: cosmoslik.plugins.Sampler
    :members:
    
