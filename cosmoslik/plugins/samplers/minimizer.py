from numpy import *
from scipy.optimize import fmin
from numpy.linalg.linalg import inv, LinAlgError
from numpy.random import multivariate_normal


def get_sampled(params): return params["_sampled"]
def get_outputted(params): return params["_output"]

def hess(f,x,dx):
    def xdx(fi,i,fj,j):
        y=list(x)
        y[abs(i)]+=fi*dx[abs(i)]
        y[abs(j)]+=fj*dx[abs(j)]
        return f(y)
    n=len(x)
    h=zeros((n,n))
    for i in range(n):
        for j in range(i,n):
            h[i,j]=h[j,i]= (xdx(1,i,1,j)-xdx(1,i,-1,j)-xdx(-1,i,1,j)+xdx(-1,i,-1,j))/(4*dx[i]*dx[j])
    return h


def sample(x,lnl,**kwargs):
    """
    Find a best fit point and possibly a Hessian 
    
    Parameters:
        start - dictionary of parameter key/values or a filename
        lnl - a function (or list of functions) of (dict:params) which returns the negative log-likelihood
        init_fn = a function (or list of functions) of (dict:params) which is called before starting the minimization
        derived_fn = a function (or list of functions) of (dict:params) which calculates and adds derived parameters 
        step_fn - a function (or list of functions) of (dict:params, dict:samples) which is called after every step
    
    Returns:
        The best fit step.
        
    """
    
    print "Minimizing..."
    def flnl(x):
        l = lnl(x,**kwargs)[0]
        print "like=%.2f step={%s}" % (l,', '.join(['%s:%.4g'%(k,v) for k,v in zip(kwargs['_sampled'],x)])) 
        return l

    if kwargs.get('minimizer.minimize',True):
        xopt, lnlopt = fmin(flnl,x,full_output=True,ftol=.01,xtol=1)[:2]
    else:
        xopt = x
        
    yield xopt, inf, 0, None
    
    if kwargs.get('minimizer.hessian',False):
        try:  
#            from numdifftools import Hessian
#            ih = Hessian(flnl,stepNom=xopt/10)(xopt)
            print "Computing hessian..."
            h = hess(flnl,xopt,[kwargs['*'+k][-1]/10 for k in kwargs['_sampled']])
            if 'minimizer.hessian_file' in kwargs:
                with open(kwargs['minimizer.hessian_file'],'w') as f:
                    f.write('# '+' '.join(get_sampled(kwargs))+'\n')
                    savetxt(f,h)
                
            try:
                ih = inv(h)
                for x in multivariate_normal(xopt,ih,size=kwargs['samples']):
                    yield x, 1, 0, None
            except LinAlgError:
                raise LinAlgError('Numerical hessian is not invertible.')
            
        except ImportError:
            print "Please install numdifftools package for Hessian."
