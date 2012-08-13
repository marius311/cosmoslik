.. _install:

==========
Installing
==========

.. warning::

    Going to figure out a better way to install, but for now...
    
CosmoSlik is written in Python, so it does not need to be compiled. You can simply
open the archive and run ``cosmoslik.py``. However, some modules are wrappers
of other codes which *do* need to be compiled. To do this, 

1. First create a file called ``Makefile.inc``, most likely by coping ``Makefile.inc.example``. 
   Then edit the file and change the various variables to ones for your platform.
   
2. Run ``cosmoslik.py --build``. This recursively runs ``make`` or ``setup.py`` on all 
   of the modules in CosmoSlik. They will be using the variables you defined in ``Makefile.inc``.
   
Requirements
============

- `Numpy <http://numpy.scipy.org/>`_
- `Scipy <http://scipy.org/>`_
- `mpi4py <http://mpi4py.scipy.org/>`_ (optional)

Troubleshooting
===============

- If a build fails, you can go directly into a module's directory and run ``make``
  or ``setup.py`` by hand.  
      
- Modules which are wrappers of Fortran code are compiled with the `F2PY <http://www.scipy.org/F2py/>`_ package. 
  The ``F2PYFLAGS`` variable in ``Makefile.inc`` specifies flags which are passed to F2PY. 
  Some useful ones include
    - ``--fcompiler=`` The name of Fortran compiler, e.g. ``gnu`` or ``intel``. Run
      ``f2py --help-fcompiler`` to list available names
    - ``--f90exec=`` A path to a Fortran compiler
    - ``--f90flags=`` Extra Fortran flags