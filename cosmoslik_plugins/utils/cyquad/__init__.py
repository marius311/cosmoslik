import os.path as osp

if osp.exists(osp.join(osp.dirname(__file__),'cyquad.so')): 
    from .cyquad import *
