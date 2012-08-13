from collections import namedtuple
import itertools
import traceback

NOMPI = False

def get_mpi():
    
    """ Return (rank,size,comm) """
    if NOMPI or len([l[2] for l in traceback.extract_stack() if l[2] == 'mpi_map']) > 1:
        (rank,size,comm) = (0,1,None)
    else:
        try:
            from mpi4py import MPI
            comm = MPI.COMM_WORLD
            (rank,size) = (comm.Get_rank(),comm.Get_size())
        except ImportError:
            (rank,size,comm) = (0,1,None)
            
    return namedtuple('mpi',['rank','size','comm'])(rank,size,comm)

def is_master(): return get_rank()==0
def get_rank(): return get_mpi().rank
def get_size(): return get_mpi().size


def mpi_consistent(value):
    """Returns the value that the root process provided."""
    comm = get_mpi().comm
    if comm==None: return value
    else: return comm.bcast(value)

def mpi_map(function,sequence,distribute=True):
    """
    A map function parallelized with MPI. If this program was called with mpiexec -n $NUM, 
    then partitions the sequence into $NUM blocks and each MPI process does the rank-th one.
    Note: If this function is called recursively, only the first call will be parallelized

    Keyword arguments:
    distribute -- If true, every process receives the answer
                  otherwise only the root process does (default=True)
    """
    (rank,size,comm) = get_mpi()
        
    sequence = mpi_consistent(sequence)

    if (size==1):
        return map(function,sequence)
    else:
        if (distribute):
            return flatten(comm.allgather(map(function, partition(sequence,size)[rank])))
        else:
            if (rank==0):
                return flatten(comm.gather(map(function, partition(sequence,size)[rank])))
            else:
                comm.gather(map(function, partition(sequence,size)[rank]))
                return []

def flatten(l):
    """Returns a list of lists joined into one"""
    return list(itertools.chain(*l))

def partition(l, n):
    """Partition l into n nearly equal sublists"""
    division = len(l) / float(n)
    return [l[int(round(division * i)): int(round(division * (i + 1)))] for i in range(n)]
