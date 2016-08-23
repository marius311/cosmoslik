#!/usr/bin/env python

# Bootstrap setuptools installation. We require setuptools >= 3.2 because of a
# bug in earlier versions regarding C++ sources generated with Cython. See:
#    https://pypi.python.org/pypi/setuptools/3.6#id171
try:
    import pkg_resources
    pkg_resources.require("setuptools >= 3.2")
except pkg_resources.ResolutionError:
    from ez_setup import use_setuptools
    use_setuptools()

import os
import errno
import fnmatch
import sys
import shlex
from distutils.sysconfig import get_config_var, get_config_vars
from subprocess import check_output, CalledProcessError, check_call
from setuptools import setup, find_packages
from setuptools.dist import Distribution
from setuptools.command.test import test as TestCommand
from distutils.command.build_clib import build_clib
from distutils.errors import DistutilsExecError
from distutils.dir_util import mkpath
from distutils.file_util import copy_file
from distutils import log

# Apple switched default C++ standard libraries (from gcc's libstdc++ to
# clang's libc++), but some pre-packaged Python environments such as Anaconda
# are built against the old C++ standard library. Luckily, we don't have to
# actually detect which C++ standard library was used to build the Python
# interpreter. We just have to propagate MACOSX_DEPLOYMENT_TARGET from the
# configuration variables to the environment.
#
# This workaround fixes <https://github.com/healpy/healpy/issues/151>.
if get_config_var('MACOSX_DEPLOYMENT_TARGET') and not 'MACOSX_DEPLOYMENT_TARGET' in os.environ:
    os.environ['MACOSX_DEPLOYMENT_TARGET'] = get_config_var('MACOSX_DEPLOYMENT_TARGET')


# If the Cython-generated C++ files are absent, then fetch and install Cython
# as an egg. If the Cython-generated files are present, then only use Cython if
# a sufficiently new version of Cython is already present on the system.
cython_require = 'Cython >= 0.16'
log.info('checking if Cython-generated files have been built')
try:
    open('healpy/src/_query_disc.cpp')
except IOError:
    log.info('Cython-generated files are absent; installing Cython locally')
    Distribution().fetch_build_eggs(cython_require)
else:
    log.info('Cython-generated files are present')
try:
    log.info('Checking for %s', cython_require)
    pkg_resources.require(cython_require)
except pkg_resources.ResolutionError:
    log.info('%s is not installed; not using Cython')
    from setuptools.command.build_ext import build_ext
    from setuptools import Extension
else:
    log.info('%s is installed; using Cython')
    from Cython.Distutils import build_ext, Extension


class custom_build_ext(build_ext):
    def finalize_options(self):
        build_ext.finalize_options(self)

        # Make sure that Numpy is importable
        # (does the same thing as setup_requires=['numpy'])
        self.distribution.fetch_build_eggs('numpy')
        # Prevent numpy from thinking it is still in its setup process:
        # See http://stackoverflow.com/questions/19919905
        __builtins__.__NUMPY_SETUP__ = False

        # Add Numpy header search path path
        import numpy
        self.include_dirs.append(numpy.get_include())

    def run(self):
        # If we were asked to build any C/C++ libraries, add the directory
        # where we built them to the include path. (It's already on the library
        # path.)
        if self.distribution.has_c_libraries():
            self.run_command('build_clib')
            build_clib = self.get_finalized_command('build_clib')
            for key, value in build_clib.build_args.items():
                for ext in self.extensions:
                    if not hasattr(ext, key) or getattr(ext, key) is None:
                        setattr(ext, key, value)
                    else:
                        getattr(ext, key).extend(value)
        build_ext.run(self)



setup(name='CosmoSlik',
      version='1.0.0',
      description='Cosmology Sampler of Likelihoods',
      classifiers=[
        #   'Development Status :: 5 - Production/Stable',
          'Environment :: Console',
          'Intended Audience :: Science/Research',
        #   'License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)',
          'Operating System :: POSIX',
          'Programming Language :: C++',
          'Programming Language :: Python :: 3.2',
          'Programming Language :: Python :: 3.3',
          'Programming Language :: Python :: 3.4',
          'Programming Language :: Python :: 3.5',
          'Topic :: Scientific/Engineering :: Astronomy',
          'Topic :: Scientific/Engineering :: Visualization'
      ],
      author='Marius Millea',
      author_email='mariusmillea@gmail.com',
      url='http://github.com/marius311/cosmoslik',
      packages=find_packages(),
      libraries=[
        #   ('cfitsio', {
        #    'pkg_config_name': 'cfitsio',
        #    'local_source': 'cfitsio',
        #    'supports_non_srcdir_builds': False}),
        #   ('healpix_cxx', {
        #    'pkg_config_name': 'healpix_cxx >= 3.30.0',
        #    'local_source': 'healpixsubmodule/src/cxx/autotools'})
      ],
      py_modules=[
        #   'healpy.pixelfunc','healpy.sphtfunc',
        #           'healpy.visufunc','healpy.fitsfunc',
        #           'healpy.projector','healpy.rotator',
        #           'healpy.projaxes','healpy.version'
                  ],
      cmdclass={
          'build_ext': custom_build_ext,
        #   'build_clib': build_external_clib,
      },
      ext_modules=[
          Extension("cosmoslik_plugins.models.cosmo_derived.cosmo_derived",
                    ['cosmoslik_plugins/models/cosmo_derived/cosmo_derived.pyx'],
                    language='c++'),
          Extension("cosmoslik_plugins.utils.cyquad.cyquad",
                    ['cosmoslik_plugins/utils/cyquad/cyquad.pyx'],
                    language='c++'),
      ],
      package_data = {
        #   'healpy': ['data/*.fits', 'data/totcls.dat', 'test/data/*.fits', 'test/data/*.sh']
          },
      scripts = ['bin/cosmoslik'],
      install_requires=['scipy', 'numpy'],
      license='GPLv2'
  )
