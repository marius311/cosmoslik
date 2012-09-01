from cosmoslik.plugins import Likelihood
from numpy import inf

class priors(Likelihood):
    def lnl(self, p, model):
        for n, lower, upper in p.get('priors',{}).get('hard',[]):
            if not (lower<p[n]<upper): return inf
        
        return sum((p[n]-c)**2/2/w**2 for n,c,w in p.get('priors',{}).get('gaussian',[]))
        