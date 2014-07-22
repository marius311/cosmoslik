ctypedef double (*cyquadfunc)(double,...)

cdef double cyquad(cyquadfunc func, double a, double b, double epsrel, double *args, int nargs)

