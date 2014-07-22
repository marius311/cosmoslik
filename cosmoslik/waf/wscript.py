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
        

def recurse(ctx, keep_going=False, exclude=[]):
    success, fail = [],[]
    for f in ctx.path.ant_glob('cosmoslik_plugins/**/wscript'):
        if not f.parent in exclude:
            try:
                ctx.recurse(f.parent.abspath())
                success.append(f.parent)
            except: 
                if keep_going: fail.append(f.parent)
                else: raise
    return success,fail

    
def options(opt, keep_going=True):
    opt.add_option('--inplace', action='store_true', default=False,
                help='install CosmoSlik in this directory')
    opt.add_option('--plugin', nargs=1, help='configure single plugin', dest='PLUGIN')
    opt.add_option('--exclude', action='append', nargs=1, help="don't build certain plugins", dest='EXCLUDE')
    recurse(opt, keep_going=True)


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


@conf
def check_lapack_blas(conf): 
    conf.env.table.setdefault('LINKFLAGS_LAPACK',['-llapack','-lblas'])
    conf.check_library_func('lapack','dpotrf_','LAPACK')
    conf.check_library_func('blas','ddot_','LAPACK')

@conf
def check_cfitsio(conf): 
    conf.env.table.setdefault('LINKFLAGS_CFITSIO',['-lcfitsio'])
    conf.check_library_func('cfitsio','ftopen_','CFITSIO')


def configure(conf):
    conf.env.PYTHON = sys.executable
    conf.load('python compiler_c')
    conf.check_python_version((2,7))
    conf.check_python_module('numpy','ver >= num(1,5)')

    for k,v in conf.environ.items():
        if any(k.startswith(p) for p in ['LIB','LIBPATH','LINKFLAGS','CFLAGS','FCFLAGS']): 
            conf.env.append_value(k, to_list(conf.environ[k]))

    if conf.options.PLUGIN is not None:
        success, fail = [conf.srcnode.find_node("cosmoslik_plugins").make_node(conf.options.PLUGIN)], []
        conf.recurse(success[0].abspath())
    else:
        if conf.options.EXCLUDE is not None:
            exclude = [conf.srcnode.find_node("cosmoslik_plugins/").find_node(e) for e in conf.options.EXCLUDE]
        else:
            exclude = []
        success, fail = recurse(conf,keep_going=True,exclude=exclude)
        
    if len(success)>0:
        sys.stdout.write('\033[1m')
        print "The following plugins are ready to build:"
        sys.stdout.write('\033[92m')
        for f in success: print "  "+f.path_from(conf.srcnode.find_node("cosmoslik_plugins/"))
        sys.stdout.write('\033[0m')
    if len(fail)>0:
        sys.stdout.write('\033[1m')
        print "The following plugins can't be built (ignore if not needed):"
        sys.stdout.write('\033[93m')
        for f in fail: print "  "+f.path_from(conf.srcnode.find_node("cosmoslik_plugins/"))
        sys.stdout.write('\033[0m')
        print "Run './waf configure --plugin PLUGIN' to see why a given plugin can't build."
        print "where PLUGIN is exactly as it appears above."

    conf.env.configured_plugins = [f.path_from(conf.srcnode) for f in success]

    if conf.options.inplace: 
        conf.env.PYTHONDIR = '.'
        conf.env.INPLACE = True
    
def build(bld):
    for f in bld.env.configured_plugins: bld.recurse(f)
    if not bld.env.INPLACE:
        bld(features='py',source=bld.path.ant_glob('cosmoslik*/**/*.py',excl='**/*waf*'),install_from='.')
    
    
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
