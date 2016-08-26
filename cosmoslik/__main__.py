#!/usr/bin/env python

import sys, os, traceback, argparse
from .cosmoslik import load_script, Slik

parser = argparse.ArgumentParser(prog='cosmoslik')
parser.add_argument('-n',type=int,default=False,help='run multiple chains with MPI')
parser.add_argument('--traceback',action='store_true',default=False,help='print out tracebacks on error messages')
parser.add_argument('script',metavar="script.py",nargs='?',help='a script to run a chain from')
parser.add_argument('script_args', nargs=argparse.REMAINDER, metavar="...", help='arguments to pass to script')

def main():

    if args.n:
        try:
            from mpi4py import MPI
        except Exception as e:
            raise Exception("Failed to load mpi4py which is needed to run with CosmoSlik '-n'.") from e
        else:
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
