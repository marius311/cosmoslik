#!/usr/bin/env python

from setuptools import setup, find_packages
from setuptools.command.sdist import sdist as _sdist
from distutils.command.build import build as _build
from distutils import file_util
import os

try:
    from sphinx.setup_command import BuildDoc as _BuildDoc
    has_sphinx = True
except ImportError:
    has_sphinx = False


class sdist(_sdist):
    """
    Modified sdist which
    * first builds the sphinx documentation so we can include it
    * follows the symlink on local_camb4py.
    """
    def run(self):
        if has_sphinx: self.run_command("build_sphinx") 
        _sdist.run(self)

    def make_release_tree(self, base_dir, files):
        _sdist.make_release_tree(self, base_dir, files)
        src_file = 'cosmoslik/models/camb/local_camb4py.py'
        target_file = os.path.join(base_dir,src_file)
        if src_file in files:
            try: os.unlink(target_file)
            except: pass
            try: self.copy_file(src_file, target_file,link=None)
            except: pass    
                    


if has_sphinx:
    class BuildDoc(_BuildDoc):
        """
        Modified build_sphinx to build the documenation in the same 
        place "make html" would, under doc/_build.
        """
        def finalize_options(self):
            self.build_dir = "doc/_build"
            _BuildDoc.finalize_options(self)
            
        def run(self):
            _BuildDoc.run(self)
            os.symlink("doc/_build/html/index.html", "README.html")


class build(_build):
    """
    Modified build to run cosmoslik.py --build, then copy the built files 
    (which get built in-place in cosmoslik/) over to build/ so that install 
    automatically sees them and installs them. 
    """
    def run(self): 
#        if has_sphinx: self.run_command("build_sphinx") 
        _build.run(self)
        #run cosmoslik --build
        #copy files to build/ directory for optional installation


cmdclass = {}
cmdclass['sdist']=sdist
cmdclass['build']=build
if has_sphinx: cmdclass['build_sphinx']=BuildDoc


setup(
    name='cosmoslik',
    version='0.1.0',
    author='Marius Millea',
    author_email='mmillea@ucdavis.edu',
    packages=find_packages(),
#    namespace_packages = ['cosmoslik','cosmoslik.plugins'],
    url='http://pypi.python.org/pypi/cosmoslik/',
    license='LICENSE.txt',
    description='A modular cosmology likelihood sampler.',
    long_description=open('README.rst').read(),
    cmdclass=cmdclass
)
