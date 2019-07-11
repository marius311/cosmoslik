#!/usr/bin/env python

import sys, os, traceback, argparse
from .cosmoslik import load_script, Slik

parser = argparse.ArgumentParser(prog='cosmoslik')
parser.add_argument('-n',type=int,default=1,help='run multiple chains with MPI')
parser.add_argument('-qn',type=int,help='submit as qsub job on QN nodes. same as --qsub "-l nodes=QN"')
parser.add_argument('--qsub',type=str,default='',help="submit as qsub job, with QSUB as options to qsub")
parser.add_argument('--qlog',type=str,default='/dev/null',help="qsub output filename (default: none)")
parser.add_argument('--qname',type=str,default='/dev/null',help="qsub job name")
parser.add_argument('--traceback',action='store_true',default=True,help='print out tracebacks on error messages')
parser.add_argument('script',metavar="script.py",help='a script to run a chain from')
parser.add_argument('script_args', nargs=argparse.REMAINDER, metavar="...", help='arguments to pass to script')

def main(args):
    
    args.qsub = args.qsub.split()
    if args.qn:
        args.qsub += ['-l','nodes=%i'%args.qn]

    if args.qsub or args.n:
        i = sys.argv.index(args.script)
        cmd = "{python} -m cosmoslik {args}".format(python=sys.executable, args=' '.join(map(escape_string,sys.argv[i:])))

    if args.qsub:
        args.qlog = args.qlog.format(**vars(args)) # allows using e.g. --qlog logs/{qname}
        args.qsub += ['-j','oe','-o',args.qlog]
        if args.qname: args.qsub += ['-N',args.qname]
        args.qsub = ' '.join(args.qsub)
        qsub_command = "echo 'cd {curdir} && mpiexec -n {n} {cmd} 2>&1 | tee {curdir}/{qlog} ' | qsub {qsub}".format(
            curdir=os.path.abspath(os.curdir), cmd=cmd, **vars(args)
        )
        print(qsub_command)
        os.system(qsub_command)
    elif args.n>1:
        try:
            from mpi4py import MPI
        except Exception as e:
            raise Exception("Failed to load mpi4py which is needed to run with CosmoSlik '-n'.") from e
        else:
            os.system("mpiexec -n {n} {cmd}".format(n=args.n, cmd=cmd))
    else:
        parser, script = load_script(args.script)
        for _ in Slik(script(**vars(parser.parse_args(args.script_args)))).sample(): pass

def escape_string(s):
    """
    Returns string with appropriate characters escaped so that it can be
    passed as a shell argument.
    """
    try:
        from subprocess import check_output
        return check_output(["bash","-c",'printf "%q" "$@"','_', s]).decode()
    except Exception as e:
        print("WARNING: `cosmoslik -n ...` may not work properly because of: \n"+str(e))
        return s

if not sys.argv[1:] or sys.argv[1] in ["-h","--help"]: 
    parser.print_help()
else:
    args = parser.parse_args()
    try:
        main(args)
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
