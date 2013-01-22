.. _install:

============
Installation
============

CosmoSlik uses the `waf <http://code.google.com/p/waf/>`_ build system. 
Installation consists of three steps::

    ./waf configure [--prefix=/my/prefix]
    ./waf build
    ./waf install
    
``configure`` verifies that all dependencies are met and searches for compilers and libraries. 
When things are not automatically configured, you can set them manually with environment variables, 
for example set the C compiler via the ``CC`` variable::

    CC=/path/to/gcc ./waf configure
    
When ``configure`` is unable to find a required
library, it will tell you the environment variables which were searched for the library. For example,
when linking to the Lapack library, one of the variables is ``LINKFLAGS_LAPACK``, which can then be
manually specified via::

    LINKFLAGS_LAPACK="-L/path/to/lapack -llapack" ./waf configure

The optional ``--prefix`` argument specifies where to install CosmoSlik after it is built. 
CosmoSlik is Python package and is installed entirely in the ``site-packages`` directory 
of a Python installation. Generally, there are two options for ``prefix``:

* left unspecified - This installs CosmoSlik system-wide, (e.g. in `/usr/lib/python/site-packages`).
  Note you may need super-user privelages for this.
* ``--prefix=~/.local`` - This installs CosmoSlik in the current user's home directory, 
  (e.g. in `~/.local/lib/python/site-packages`)

Once the ``configure`` step finishes successfully, ``build`` is ready to run. This creates a folder
``build/`` which houses the compiled libraries.

  


Requirements
============

- `Numpy <http://numpy.scipy.org/>`_
- `Scipy <http://scipy.org/>`_
- `mpi4py <http://mpi4py.scipy.org/>`_ (optional)
