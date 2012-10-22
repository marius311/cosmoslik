APPNAME = 'cosmoslik'
VERSION = '0.1'


#def dist(ctx):
#    ctx.algo = 'tar.gz' 
#    #ctx.excl = ' **/.waf-1* **/*~ **/*.pyc **/*.swp **/.lock-w* **/.git*' 

def recurse(ctx):
    for f in ctx.path.ant_glob('cosmoslik/**/wscript'): ctx.recurse(f.parent.abspath())
    
def options(opt):
    recurse(opt)
    
def configure(conf):
    recurse(conf)

    
def build(bld):
    recurse(bld)
    pyfiles = [f.relpath() for f in bld.path.ant_glob('cosmoslik/**/*.py',excl='**/*waf*')]
    bld(features='subst',source=pyfiles,target=pyfiles,is_copy=True)



