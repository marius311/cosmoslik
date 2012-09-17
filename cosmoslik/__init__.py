#The magic line to make this a "namespace" package.
#see http://packages.python.org/distribute/setuptools.html#namespace-packages
__path__ = __import__('pkgutil').extend_path(__path__,__name__)

from cosmoslik import *
import params, chains
