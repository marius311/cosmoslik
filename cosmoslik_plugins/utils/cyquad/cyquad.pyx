#cython: profile=False

import cython
from libc.math cimport INFINITY as inf


#steal pointers directly to scipy's DQAGSE and DQAGIE fortran functions
import ctypes, scipy.integrate

qp=ctypes.cdll.LoadLibrary(scipy.integrate._quadpack.__file__)

ctypedef void (*dqagse_fn)(double (*f)(double *), double *a, double *b, double *epsabs, double *epsrel, 
                           int *limit, double *result, double *abserr, int *neval, int *ier,
                           double *alist, double *blist, double *rlist, double *elist, int *iord, int *last)
cdef dqagse_fn dqagse = (<dqagse_fn*><size_t>ctypes.addressof(qp.dqagse_))[0]

ctypedef void (*dqagie_fn)(double (*f)(double *), double *bound, int *inf, double *epsabs, double *epsrel, 
                           int *limit, double *result, double *abserr, int *neval, int *ier, 
                           double *alist, double *blist, double *rlist, double *elist, int *iord, int *last)
cdef dqagie_fn dqagie = (<dqagie_fn*><size_t>ctypes.addressof(qp.dqagie_))[0]
#----



#keep a call stack to pass *args to the integrand function, since DQAGSE doesn't allow this directly
cdef struct call_sig:
    cyquadfunc func
    double *args
    int nargs
    
cdef call_sig init_call_sig(cyquadfunc func, double *args, int nargs):
    cdef call_sig sig
    sig.func=func
    sig.args=args
    sig.nargs=nargs
    return sig

cdef class call_stack:
    #hats off to you if you break my crappy stack implementation by doing >100 dimensional integration
    cdef call_sig sigs[100] 
    cdef int top

    cdef void push(self, call_sig sig):
        self.sigs[self.top] = sig
        self.top = self.top + 1
    
    cdef call_sig pop(self):
        self.top = self.top-1
        return self.sigs[self.top]        
    
    cdef inline call_sig peek(self):
        return self.sigs[self.top-1]        

cdef call_stack my_call_stack = call_stack()
#---- 

@cython.profile(False)
cdef double quad_func_wrapper(double *x):
    cdef call_sig sig = my_call_stack.peek()
    cdef double *a = sig.args

    #ugly and there's a max argument limit, but hey if this part's hidden the rest is beautiful....   
    if sig.nargs==0:
        return sig.func(x[0])
    elif sig.nargs==1:
        return sig.func(x[0],a[0])
    elif sig.nargs==2:
        return sig.func(x[0],a[0],a[1])
    elif sig.nargs==3:
        return sig.func(x[0],a[0],a[1],a[2])
    elif sig.nargs==4:
        return sig.func(x[0],a[0],a[1],a[2],a[3])
    elif sig.nargs==5:
        return sig.func(x[0],a[0],a[1],a[2],a[3],a[4])
    elif sig.nargs==6:
        return sig.func(x[0],a[0],a[1],a[2],a[3],a[4],a[5])
    elif sig.nargs==7:
        return sig.func(x[0],a[0],a[1],a[2],a[3],a[4],a[5],a[6])
    elif sig.nargs==8:
        return sig.func(x[0],a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7])
    elif sig.nargs==9:
        return sig.func(x[0],a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8])
    elif sig.nargs==10:
        return sig.func(x[0],a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9])





cdef double cyquad(cyquadfunc func, double a, double b, double epsrel, double *args, int nargs):
    """
    Call SciPy's quad function (DQAGSE) directly from Cython, with no Python function call overhead. 
    
    Args:
        cyquadfunc func : pointer to the integrand function
        double a,b : limits of integration (can be inf's)
        double *args : an array of extra arguments to pass to the integration function (max length: 10)
        int nargs : the length of args
        
    Example:

        >> cdef double f(double x, double m, double b):
        >>    return m*x+b
        >> cyquad(<cyquadfunc>f,0,1,[3,1],2) #integrates 3*x+1 from 0 to 1
        2.5

    Notes:
        * multiple integration is possible, so `func` can itself call cyquad recursively
        * this implementation is NOT thread safe

    """
    
    cdef double epsabs=0, bound=0, result, abserr
    cdef int limit=50, infbound=1, neval,ier,last  
    cdef double alist[50], blist[50], rlist[50], elist[50]
    cdef int iord[50]

    my_call_stack.push(init_call_sig(func,args,nargs))
    if a!=-inf and b!=inf:
        dqagse(quad_func_wrapper,&a,&b,&epsabs,&epsrel,&limit,&result,&abserr,&neval,
               &ier,alist,blist,rlist,elist,iord,&last)
    else:
        if a==-inf and b==inf: infbound, bound = 2, 0
        elif a==-inf: infbound, bound = -1, b
        else: infbound, bound = 1, a
        dqagie(quad_func_wrapper,&bound,&infbound,&epsabs,&epsrel,&limit,&result,&abserr,&neval,
           &ier,alist,blist,rlist,elist,iord,&last)
    my_call_stack.pop()
    # print result, (abserr/result)/epsrel, neval, ier, [args[i] for i in range(nargs)], args[1]/args[0]

    return result
