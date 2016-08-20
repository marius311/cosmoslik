#!/usr/bin/env python

import sys, os, traceback, argparse
from .cosmoslik import get_plugin, get_all_plugins, load_script, Slik

parser = argparse.ArgumentParser(prog='cosmoslik')
parser.add_argument('-n',type=int,default=False,help='run multiple chains with MPI')
parser.add_argument('--traceback',action='store_true',default=False,help='print out tracebacks on error messages')
group = parser.add_mutually_exclusive_group()
group.add_argument('--list',action='store_true',default=False,help='list available modules')
group.add_argument('--doc',nargs=1,metavar='<module>',help='print the documentation for a module')
group.add_argument('--html_doc',nargs=1,metavar='<module>',help='open the documentation for a module in a web-browser')
group.add_argument('script',metavar="script.py",nargs='?',help='a script to run a chain from')
parser.add_argument('script_args', nargs=argparse.REMAINDER, metavar="...", help='arguments to pass to script')

def main():

    if args.list:
        print("Found the following modules in 'cosmoslik_plugins':")
        for name in sorted(get_all_plugins().values()):
            print('  %s'%'.'.join(name.split('.')[1:]))
        print("See 'cosmoslik.py --doc <module>' for more information on a given module.")
        print("Some modules may need to be compiled before appearing in this list.")
        
    elif args.doc or args.html_doc:
        from textwrap import dedent
        plugin_name = (args.doc or args.html_doc)[0]
        plugin = get_plugin(plugin_name) 
        doc = plugin.__doc__ or ""
        if args.doc:
            print("Documentation for module '%s':"%plugin_name)
            print(dedent(doc))
        else:
            from docutils.core import publish_string
            from tempfile import mktemp
            import webbrowser
            tmpfile = mktemp(suffix='.html')
            with open(tmpfile,'w') as f: f.write(publish_string(dedent(doc),writer_name='html'))
            webbrowser.open(tmpfile)

    elif args.n:
        i = sys.argv.index('-n')
        sys.argv.pop(i); sys.argv.pop(i)
        os.system("mpiexec -n %i %s -m cosmoslik %s"%(args.n,sys.executable,' '.join(sys.argv[1:])))
        
    elif args.script:
        parser, script = load_script(args.script)
        for _ in Slik(script(**vars(parser.parse_args(args.script_args)))).sample(): pass



if not sys.argv[1:]: parser.print_help()
else:
    args = parser.parse_args()
    try:
        main()
    except BaseException as e:
        if isinstance(e,(Exception,KeyboardInterrupt)):
            if args.traceback: 
                traceback.print_exception(type(e), e, sys.exc_info()[2], None, sys.stderr)
            else:
                sys.stderr.write('\033[91m')
                traceback.print_exception(type(e), e, None, None, sys.stderr)
                sys.stderr.write('\033[0m')
                print("Run CosmoSlik with --traceback for more info.")
            sys.exit(1)
        else:
            raise

    
    
