#!/usr/bin/env python

import sys, os, cosmoslik, traceback, argparse

parser = argparse.ArgumentParser(prog='cosmoslik.py')
parser.add_argument('params.ini',nargs='*',help='a parameter file to run')
parser.add_argument('--list',action='store_true',default=False,help='list available modules')
parser.add_argument('--doc',nargs=1,metavar='<module>',help='print the documentation for a module')
parser.add_argument('--html_doc',nargs=1,metavar='<module>',help='open the documentation for a module in a web-browser')
parser.add_argument('--build',nargs='?',metavar='<modules>',default=False,help='run build script for a module (default: all modules)')
parser.add_argument('--clean',nargs='?',metavar='<modules>',default=False,help='run clean for a module (default: all modules)')
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
        
    elif args['clean'] is not False:
        cosmoslik.build(args['clean'],'clean',('cleaned','clean'))

    elif args['qsub']:
        from subprocess import Popen, PIPE
        from cosmoslik import params

        inifiles = args['params.ini']
        
        if ('camb' in params.read_ini(inifiles[0]).get('models',[])): nnodes, ppn, wall, nproc = 6, 1, 24, 7
        else: nnodes, ppn, wall, nproc = 1, 1, 6, 9
        #TODO: allow mixing CAMB/PICO runs
        nnodes *= len(inifiles)
            
        name = 'cosmoslik'
        sys.argv.remove('--qsub')
        
        proc = Popen(["qsub","-q","usplanck",
                      "-l","nodes=%s:ppn=%s,pvmem=20gb"%(nnodes,ppn),
                      "-l","walltime=%s:00:00"%wall,
                      "-N",name,"-o","%s.log"%name,"-j","oe","-V"],stdin=PIPE,stdout=PIPE)
        
        cmd = ' &\n'.join(['(cd %s && %s -m cosmoslik -n %i %s &> %s)'%\
                           (os.path.dirname(os.path.abspath(inifile)),
                            sys.executable,
                            nproc,
                            os.path.basename(inifile),
                            os.path.basename(inifile).replace('.ini','')+'.log') for inifile in inifiles])
        print cmd
        proc.stdin.write(cmd)
        proc.stdin.close()
        print proc.stdout.readline()

    elif args['n']:
        i = sys.argv.index('-n')
        sys.argv.pop(i); sys.argv.pop(i)
        os.system("mpiexec -n %i %s -m cosmoslik %s"%(int(args['n'][0]),sys.executable,' '.join(sys.argv[1:])))
        
    elif args['params.ini']:
        if len(args['params.ini'])>1:
            raise Exception("Running multiple parameter files at once only available with --qsub.")
        else:
            for _ in cosmoslik.sample(args['params.ini'][0]): pass



if not sys.argv[1:]: parser.print_help()
else:
    
    args = vars(parser.parse_args())
    try:
        main(args)
    except BaseException as e:
        sys.stderr.write('\033[91m')
        traceback.print_exception(type(e), e, None, None, sys.stderr)
        sys.stderr.write('\033[0m')
        if args['traceback']: traceback.print_exception(None, None, sys.exc_info()[2], None, sys.stderr)
        else: print "Run CosmoSlik with --traceback for more info."
        sys.exit(1)

    
    
