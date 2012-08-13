============
The Ini File
============

Format
======

- CosmoSlik ini files are lists of ``key = value`` pairs.
- White space and indentation are ignored.
- They can have *sections*
- The comment character is #
- Values which are valid Python expression are dynamically evaluated

For example::

    file = /some/path
    
    #A comment
    frequencies = ['143','217']
     


Sections
--------

For organizational purposes, CosmoSlik supports sections. The syntax is similar to 
Windows ini files, but can be recursive. ::

    [section]{
        key1 = value
        
        [sub_section]{
            # a sampled parameter in a subsection
            key2 = 0 [-10 10 1]
        }
    }

When parameters are sampled in a subsection, their name appears in the output file 
joined by a ``.`` (period). For example, here it would be called ``section.sub_section.key2``. 

Sampling Parameters
-------------------

To have the sampler vary and sample a parameter, add ``[min max width]`` after its value. 
For example::

    H0 = 70 [40 100 5]
    
will sample over a parameter named ``H0`` restricted to the range ``40<H0<70`` and with 
a natural width of ``5`` (different samplers may use this value differently).



Including other Files
---------------------

CosmoSlik supports *including* other ini files via::

    include other_file.ini
    
It will be exactly as if the contents of the other file are pasted at the location
of the ``inculde``. 

Recursive includes are not supported.


Parameters
==========

These parameters are built in to CosmoSlik. Each module adds other parameters. 
(see individual documentation). 

update_frequency
----------------
    Print chain info to the screen every # of steps. (default: 1)
    
output_file
-----------
    The chain gets outputted to this file 

samples
-------
    The number of (possibly non-unique) samples to draw. (default: 100)
    
    
likelihooods
------------
models
------
samplers
--------
derivers
--------
    Each of these keys should be given a list of names of modules to dynamically load.
    For available modules, see documentation, or run ``cosmoslik.py --list``. For example,
    if ::
    
        $ ./cosmoslik.py --list
        Found the following modules:
            likelihoods.some_module
            
    then ::
    
        likelihoods = some_module
        
    would be a valid line in the ini file.
        
