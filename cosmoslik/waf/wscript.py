#! /usr/bin/env python

APPNAME = 'cosmoslik'
VERSION = '0.1'


import sys, os.path as osp, os
from waflib.Utils import to_list
from waflib.Configure import conf
from ConfigParser import RawConfigParser
from cStringIO import StringIO
from setuptools.command import easy_install

#def dist(ctx):
#    ctx.algo = 'tar.gz' 
#    ctx.excl = ' **/.waf-1* **/*~ **/*.pyc **/*.swp **/.lock-w* **/.git*' 

#def update_enivorn(ctx,ini):
#    ini = open(ini).read()
#    config = RawConfigParser()
#    config.optionxform=str
#    config.readfp(StringIO('[root]\n'+ini))
#    for k,v in config.items('root'):
#        ctx.environ[k.upper()] = v
        

def recurse(ctx):
    for f in ctx.path.ant_glob('cosmoslik_plugins/**/wscript'): 
        ctx.recurse(f.parent.abspath())
    
def options(opt):
    recurse(opt)
    

@conf
def check_library_func(conf, library, function, use, envvars=None):
    if envvars is None: envvars = []
    envvars += ['LIB','LIBPATH','LINKFLAGS']
    for u in to_list(use): envvars += ['LIB_%s'%u,'LIBPATH_%s'%u,'LINKFLAGS_%s'%u]
    try:
        conf.check(fragment='int main(){ %s(); }'%function,
                   msg="Checking for library %s"%library,
                   okmsg="found function %s"%function,
                   errmsg="couldn't find function %s"%function,
                   use=use)
    except conf.errors.ConfigurationError: 
        conf.fatal(('%s was not found on your system. '
                    'Check that it is installed and that the following '
                    'environment variables are set correctly:\n')%library+
                    '\n'.join(['%s = %s'%(x,' '.join(getattr(conf.env,x,''))) for x in sorted(set(envvars))]))



def configure(conf):
    conf.env.PYTHON = sys.executable
#    if osp.exists('setup.ini'): update_enivorn(conf, 'setup.ini')
    
    conf.load('python compiler_c')
    conf.check_python_version((2,7))
    conf.check_python_module('numpy','ver >= num(1,5)')
    try:
        conf.check_python_module('pypico','ver == num(3,1,0)')
    except: 
        args =sys.argv[2:]+["-U","pypico==3.1.0"]
        print args
        easy_install.main(args)

    
    for x,d in [('LINKFLAGS',None),
                ('LINKFLAGS_LAPACK','-llapack -lblas'),
                ('LINKFLAGS_CFITSIO','-lcfitsio')]:
        conf.env.append_value(x, to_list(conf.environ.get(x,d or [])))

    for lib, func, use in [('lapack','dpotrf_','LAPACK'),
                           ('blas','ddot_','LAPACK'),
                           ('cfitsio','ftopen_','CFITSIO')]:
        conf.check_library_func(lib,func,use=use)

    recurse(conf)

    
def build(bld):
    recurse(bld)
    bld(features='py',source=bld.path.ant_glob('cosmoslik*/**/*.py',excl='**/*waf*'),install_from='.')



def develop(bld):
    from setuptools import setup, find_packages
    setup(
        name=APPNAME,
        version=VERSION,
        author='Marius Millea',
        author_email='mmillea@ucdavis.edu',
        packages=find_packages(),
        description='A cosmologial sampler of likelihooods.',
    )

    
"""
Subclass C/FC program so that LINKFLAGS get put at the end
"""
import waflib.Tools.c, waflib.Tools.fc 
class cprogram(waflib.Tools.c.cprogram):
    run_str='${LINK_CC} ${CCLNK_SRC_F}${SRC} ${CCLNK_TGT_F}${TGT[0].abspath()} ${RPATH_ST:RPATH} ${FRAMEWORKPATH_ST:FRAMEWORKPATH} ${FRAMEWORK_ST:FRAMEWORK} ${ARCH_ST:ARCH} ${STLIB_MARKER} ${STLIBPATH_ST:STLIBPATH} ${STLIB_ST:STLIB} ${SHLIB_MARKER} ${LIBPATH_ST:LIBPATH} ${LIB_ST:LIB} ${LINKFLAGS}'
    
class cshlib(cprogram,waflib.Tools.c.cshlib): pass

class fcprogram(waflib.Tools.fc.fcprogram):
    run_str='${FC} ${FCLNK_SRC_F}${SRC} ${FCLNK_TGT_F}${TGT[0].abspath()} ${RPATH_ST:RPATH} ${FCSTLIB_MARKER} ${FCSTLIBPATH_ST:STLIBPATH} ${FCSTLIB_ST:STLIB} ${FCSHLIB_MARKER} ${FCLIBPATH_ST:LIBPATH} ${FCLIB_ST:LIB} ${LINKFLAGS}'
    
class fcshlib(fcprogram,waflib.Tools.fc.fcshlib): pass
