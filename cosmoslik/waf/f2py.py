import os.path as osp, numpy
from waflib.Utils import to_list
from utils import to_nodes


def build_f2py(bld, source, module_name, extra_sources, skip=None, only=None, symlink=False, **kwargs):
    """
    Build an f2py extension with waf.
    
    Arguments:
    ----------
    
    source - the name of the file being wrapped
    module_name - the name of the module being produced
    extra-sources - other things to compile and link along with the extension
    skip/only - skip/only wrap certain functions
    symlink - symlink the library into the source folder after building
    **kwargs - passed to the build command
    """
    
    #use f2py to create the wrapper
    skip = 'skip: %s'%' '.join(to_list(skip)) if skip is not None else ''
    only = 'only: %s'%' '.join(to_list(only)) if only is not None else ''
    bld(rule=('${{F2PY}} --build-dir ${{TGT[0].parent.abspath()}} --quiet '
              '-m {MODULENAME} ${{SRC}} {ONLY} {SKIP}').format(MODULENAME=module_name,
                                                               ONLY=only,
                                                               SKIP=skip), 
        source=source, 
        target='{MODULENAME}module.c {MODULENAME}-f2pywrappers2.f90'.format(MODULENAME=module_name))

    #copy over f2py files fortranobject.c and fortranobject.h
    import numpy.f2py.cfuncs as cfuncs
    f2pydir = osp.join(osp.dirname(osp.abspath(cfuncs.__file__)),'src')
    for x in 'ch':
        bld(rule='cp -L ${SRC} ${TGT}',
            source=bld.root.find_node(osp.join(f2pydir,'fortranobject.'+x)),
            target='fortranobject.'+x)
    
    #build the library
    bld(features='c fc fcshlib', 
        source=('{SRC} '
                'fortranobject.c '
                '{MODULENAME}module.c '
                '{MODULENAME}-f2pywrappers2.f90').format(SRC=source,
                                                         MODULENAME=module_name).split() + to_nodes(bld,extra_sources),
        target=module_name, 
        includes=[bld.root.find_node(numpy.get_include())]+to_nodes(bld,kwargs.pop('includes',[])),
        use=['PYEXT']+to_list(kwargs.pop('use',[])),
        install_path=kwargs.pop('install_path',None) or osp.join('${PYTHONDIR}',bld.path.path_from(bld.root.find_node(bld.top_dir))),
        **kwargs)
    
    if symlink:
        bld.symlink_as(bld.path.get_src().make_node('%s.so'%module_name).abspath(),
                       bld.path.get_bld().make_node('%s.so'%module_name).abspath())
    

def config_f2py(conf):
    conf.load('compiler_c compiler_fc python')
    conf.check(features='c', cflags=['-fPIC'])
    conf.check(features='fc', cflags=['-fPIC'])
    conf.env.CFLAGS = ['-fPIC']
    conf.env.FCFLAGS = ['-fPIC']
    conf.check_python_headers()
    conf.env.fcshlib_PATTERN = '%s.so'
    conf.find_program('f2py',var='F2PY')
    
def opt_f2py(opt):
    opt.load('compiler_c compiler_fc python')

config = config_f2py
opt = opt_f2py
