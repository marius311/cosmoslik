import os.path as osp
from waflib.Utils import to_list

def fpreproc(bld,source):
    """
    Pre-process some Fortran files using the -E option before feeding them into waf.
    """
    ppfs = []
    for f in to_list(source): 
        if isinstance(f,str): f=bld.path.make_node(f)
        if f.suffix()=='.F90': 
            ppf = f.change_ext('.f90')
            bld(rule='${FC} -E ${SRC} > ${TGT}', source=f, target=ppf)
            ppfs.append(ppf)
        else:
            ppfs.append(f)
    return ppfs if len(ppfs)>1 else ppfs[0]


def build_f2py(bld, source, module_name, extra_sources):
    """
    Build an f2py extension with waf.
    
    Arguments:
    ----------
    
    source - the name of the file being wrapped
    module_name - the name of the module being produced
    extra-sources - other things to compile along with the extension
    """
    
    #use f2py to create the wrapper
    bld(rule='${{F2PY}} ${{SRC}} -m {MODULENAME} --build-dir ${{TGT[0].parent.abspath()}}'.format(MODULENAME=module_name), 
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
                                                         MODULENAME=module_name).split() + to_list(extra_sources),
        target=module_name, 
        use='PYEXT')

def config_f2py(conf):
    conf.load('compiler_c compiler_fc')
    conf.check(features='c', cflags=['-fPIC'])
    conf.check(features='fc', cflags=['-fPIC'])
    conf.env.CFLAGS = ['-fPIC']
    conf.env.FCFLAGS = ['-fPIC']
    conf.load('python')
    conf.check_python_headers()
    conf.env.fcshlib_PATTERN = '%s.so'
    conf.find_program('f2py',var='F2PY')
    
def opt_f2py(opt):
    opt.load('compiler_c compiler_fc')
    opt.load('python')

