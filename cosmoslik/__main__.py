#!/usr/bin/env python

import sys, os, cosmoslik, traceback, argparse

parser = argparse.ArgumentParser(prog='cosmoslik.py')
parser.add_argument('params.ini',nargs='?',help='a parameter file to run')
parser.add_argument('--list',action='store_true',default=False,help='list available modules')
parser.add_argument('--doc',nargs=1,metavar='<module>',help='print the documentation for a module')
parser.add_argument('--html_doc',nargs=1,metavar='<module>',help='open the documentation for a module in a web-browser')
parser.add_argument('--build',nargs='?',metavar='<modules>',default=False,help='run build script for a module (default: all modules)')
parser.add_argument('-n',nargs=1,metavar='<# of chains>',default=False,help='run multiple chains with MPI')
parser.add_argument('--qsub',action='store_true',default=False,help='submit via qsub')
parser.add_argument('--traceback',action='store_true',default=False,help='print out tracebacks on error messages')

def main(args):
    if args['list']:
        import plugins
        print "Found the following modules in 'cosmoslik.plugins':"
        for (name,_,typ) in plugins.get_all_plugins():
            print '  %s'%'.'.join(name.split('.')[2:])
        print "See 'cosmoslik.py --doc <module>' for more information on a given module."
        
    elif args['doc'] or args['html_doc']:
        from textwrap import dedent
        import plugins
        plugin_name = (args['doc'] or args['html_doc'])[0]
        plugin = plugins.get_plugin(plugin_name) 
        doc = plugin.__doc__ or ""
        if args['doc']:
            print "Documentation for module '%s':"%plugin_name
            print dedent(doc)
        else:
            from docutils.core import publish_string
            from tempfile import mktemp
            import webbrowser
            tmpfile = mktemp(suffix='.html')
            with open(tmpfile,'w') as f: f.write(publish_string(dedent(doc),writer_name='html'))
            webbrowser.open(tmpfile)

            
    elif args['build'] is not False:
        cosmoslik.build(args['build'])
        
    elif args['qsub']:
        from subprocess import Popen, PIPE
        from cosmoslik import params

        inifile = args['params.ini']
        
        if ('camb' in params.read_ini(inifile).get('models',[])): nnodes, ppn, wall, nproc = 6, 1, 24, 7
        else: nnodes, ppn, wall, nproc = 1, 1, 6, 9
        
        name = inifile.replace('cosmoslik.','').replace('params.','').replace('.ini','')
        dir = os.path.dirname(os.path.abspath(inifile))
        sys.argv.remove('--qsub')
        
        proc = Popen(["qsub","-q","usplanck","-l",
                      "nodes=%s:ppn=%s,pvmem=20gb"%(nnodes,ppn),
                      "-l","walltime=%s:00:00"%wall,
                      "-N",name,"-o","%s.log"%name,"-j","oe","-V"],stdin=PIPE,stdout=PIPE)
        
        proc.stdin.write('cd %s && %s -m cosmoslik -n %i %s'%(dir,sys.executable,nproc,' '.join(sys.argv[1:])))
        proc.stdin.close()
        print proc.stdout.readline()

    elif args['n']:
        i = sys.argv.index('-n')
        sys.argv.pop(i); sys.argv.pop(i)
        os.system("mpiexec -n %i %s -m cosmoslik %s"%(int(args['n'][0]),sys.executable,' '.join(sys.argv[1:])))
        
    elif args['params.ini']:
        for _ in cosmoslik.sample(args['params.ini']): pass



if not sys.argv[1:]: parser.print_help()
else:
    
    args = vars(parser.parse_args())
    try:
        main(args)
    except Exception as e:
        sys.stderr.write('\033[91m')
        traceback.print_exception(type(e), e, None, None, sys.stderr)
        sys.stderr.write('\033[0m')
        if args['traceback']: traceback.print_exception(None, None, sys.exc_info()[2], None, sys.stderr)
        else: print "Run CosmoSlik with --traceback for more info."
        sys.exit(1)

    
    
