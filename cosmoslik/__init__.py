from .cosmoslik import *
from . import utils

models      = cosmoslik.plugin_getter("cosmoslik_plugins.models")
likelihoods = cosmoslik.plugin_getter("cosmoslik_plugins.likelihoods")
samplers    = cosmoslik.plugin_getter("cosmoslik_plugins.samplers")
