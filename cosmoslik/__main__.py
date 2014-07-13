#!/usr/bin/env python

import sys, os, cosmoslik as K, traceback, argparse

parser = argparse.ArgumentParser(prog='cosmoslik.py')
parser.add_argument('params.ini',nargs='?',help='a parameter file to run')
parser.add_argument('--list',action='store_true',default=False,help='list available modules')
parser.add_argument('--doc',nargs=1,metavar='<module>',help='print the documentation for a module')
parser.add_argument('--html_doc',nargs=1,metavar='<module>',help='open the documentation for a module in a web-browser')
parser.add_argument('-n',nargs=1,metavar='<# of chains>',default=False,help='run multiple chains with MPI')
parser.add_argument('--traceback',action='store_true',default=False,help='print out tracebacks on error messages')

def main(args):

    if args['list']:
        print "Found the following modules in 'cosmoslik_plugins':"
        for name in sorted(K.get_all_plugins().values()):
            print '  %s'%'.'.join(name.split('.')[1:])
        print "See 'cosmoslik.py --doc <module>' for more information on a given module."
        
    elif args['doc'] or args['html_doc']:
        from textwrap import dedent
        plugin_name = (args['doc'] or args['html_doc'])[0]
        plugin = K.get_plugin(plugin_name) 
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

    elif args['n']:
        i = sys.argv.index('-n')
        sys.argv.pop(i); sys.argv.pop(i)
        os.system("mpiexec -n %i %s -m cosmoslik %s"%(int(args['n'][0]),sys.executable,' '.join(sys.argv[1:])))
        
    elif args['params.ini']:
        from cosmoslik import load_script
        p = load_script(args['params.ini'])
        for _ in p.sample(): pass



if not sys.argv[1:]: parser.print_help()
else:
    
    args = vars(parser.parse_args())
    try:
        main(args)
    except BaseException as e:
        if args['traceback']: 
            traceback.print_exception(type(e), e, sys.exc_info()[2], None, sys.stderr)
        else:
            sys.stderr.write('\033[91m')
            traceback.print_exception(type(e), e, None, None, sys.stderr)
            sys.stderr.write('\033[0m')
            print "Run CosmoSlik with --traceback for more info."
            
        sys.exit(1)

    
    
