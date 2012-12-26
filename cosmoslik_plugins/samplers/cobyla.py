from cosmoslik.plugins import Sampler
from scipy.optimize import *
from numpy import *
from numpy.linalg.linalg import inv

class cobyla(Sampler):
    
    def sample(self, x, lnl, p):
        
        
        scale = array([v[-1] for v in p.get_all_sampled().values()])
        
        def lnl_scaled(x):
            print dict(zip(['.'.join(k) for k in p.get_all_sampled().keys()],x*scale))
            return lnl(x*scale,p)[0]
        
        def lnl_unscaled(x):
            print dict(zip(['.'.join(k) for k in p.get_all_sampled().keys()],x))
            return lnl(x,p)[0]

        cons = []
        for (i,v) in enumerate(p.get_all_sampled().values()):
            cons.append(lambda x: x[i]*scale[i] - v[1])
            cons.append(lambda x: -x[i]*scale[i] + v[2])
        
        xopt = fmin_cobyla(lnl_scaled, 
                            x/scale, 
                            cons,
#                            full_output=True,
                            rhobeg=10,
                            rhoend=1e-2,
                            disp=1)*scale
        
        yield xopt, 1, lnl(xopt,p)[0], lnl(xopt,p)[1]
        
#        def hess(f,x,dx):
#            def xdx(fi,i,fj,j):
#                y=list(x)
#                y[abs(i)]+=fi*dx[abs(i)]
#                y[abs(j)]+=fj*dx[abs(j)]
#                return f(y)
#            n=len(x)
#            h=zeros((n,n))
#            for i in range(n):
#                for j in range(i,n):
#                    h[i,j]=h[j,i]= (xdx(1,i,1,j)-xdx(1,i,-1,j)-xdx(-1,i,1,j)+xdx(-1,i,-1,j))/(4*dx[i]*dx[j])
#            return h
#
#        H = hess(lnl_unscaled,xopt,scale/100)
#        print H
#            
##        from numdifftools import Hessian
##        
##        iH = inv(Hessian(lnl_unscaled,stepNom=scale/100.)(xopt))
#        print dict(zip(['.'.join(k) for k in p.get_all_sampled().keys()],sqrt(diag(inv(H)))))
