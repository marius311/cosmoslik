#!/usr/bin/env python

import sys, platform
assert sys.version_info >= (3,3), "CosmoSlik requires Python version >= 3.3. You have "+platform.python_version()

import pkg_resources
pkg_resources.require("setuptools >= 3.2")
from setuptools import setup, find_packages

pkg_resources.require('Cython >= 0.16')
from Cython.Distutils import build_ext, Extension
from Cython.Build import cythonize

setup(name='CosmoSlik',
      version='1.0.0',
      description='Cosmology Sampler of Likelihoods',
      classifiers=[
          'Development Status :: 4 - Beta',
          'Environment :: Console',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
          'Operating System :: POSIX',
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
      cmdclass={
          'build_ext': build_ext,
      },
      ext_modules=cythonize([
          Extension("cosmoslik_plugins.misc.cyquad.cyquad",
                    ['cosmoslik_plugins/misc/cyquad/cyquad.pyx'],
                    language='c'),
          Extension("cosmoslik_plugins.models.cosmo_derived.cosmo_derived",
                    ['cosmoslik_plugins/models/cosmo_derived/cosmo_derived.pyx'],
                    language='c'),
      ]),
      scripts = ['bin/cosmoslik'],
      install_requires=['scipy>=0.17.0', 'numpy>=1.11.0'],
      license='GPLv3'
  )
