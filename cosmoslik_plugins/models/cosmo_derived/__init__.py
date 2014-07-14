import os.path as osp

if osp.exists(osp.join(osp.dirname(__file__),'cosmo_derived.so')): 
    from cosmo_derived import *
