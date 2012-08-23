"""

Python module for dealing with MCMC chains.

Contains utilities to:

* Load chains in a variety of formats
* Compute statistics (mean, std-dev, conf intervals, ...)
* Plot 1- and 2-d distributions of one or multiple parameters.
* Post-process chains 

"""


import os, sys, re
from numpy import *
from itertools import takewhile
from cosmoslik.params import load_ini

__all__ = ['Chain','Chains',
           'like1d','like2d','likegrid',
           'get_covariance', 'load_chain']


class Chain(dict):
    """
    An MCMC chain. This is just a Python dictionary mapping parameter names
    to arrays of values, along with the special keys 'lnl' and 'weight'
    """
    def __init__(self,*args,**kwargs):
        super(Chain,self).__init__(*args,**kwargs)
        for k,v in self.items(): self[k]=atleast_1d(v)
        if self and 'weight' not in self: self['weight']=ones(len(self.values()[0]))
        
    def copy(self):
        """Deep copy the chain so post-processing, etc... works right"""
        return Chain({k:v.copy() for k,v in self.iteritems()})
        
    def params(self): 
        """Returns the parameters in this chain (i.e. the keys except 'lnl' and 'weight'"""
        return set(self.keys())-set(["lnl","weight"])
    
    def sample(self,s,keys=None): 
        """Return a sample or a range of samples depending on if s is an integer or a slice object."""
        return Chain((k,self[k][s]) for k in (keys if keys else self.keys()))
    
    def matrix(self,params=None):
        """Return this chain as an nsamp * nparams matrix."""
        return vstack([self[p] for p in (params if params else self.params())]).T
    
    def cov(self,params=None): 
        """Returns the covariance of the parameters (or some subset of them) in this chain."""
        return get_covariance(self.matrix(params), self["weight"])
    
    def mean(self,params=None): 
        """Returns the mean of the parameters (or some subset of them) in this chain."""
        return average(self.matrix(params),axis=0,weights=self["weight"])
    
    def std(self,params=None): 
        """Returns the std of the parameters (or some subset of them) in this chain."""
        return sqrt(average((self.matrix(params)-self.mean(params))**2,axis=0,weights=self["weight"]))
    
    def acceptance(self): 
        """Returns the acceptance ratio."""
        return 1./mean(self["weight"])
    
    def length(self,unique=True):
        """Returns the number of unique samples. Set unique=False to get total samples."""
        return (len if unique else sum)(self['weight'])
    
    def burnin(self,n):
        """Remove the first n non-unique samples from the beginning of the chain."""
        return self.sample(slice(sum(1 for _ in takewhile(lambda x: x<n, cumsum(self['weight']))),None))
    
    def best_fit(self):
        """Get the best fit sample."""
        return {k:v[0] for k,v in self.sample(self['lnl'].argmin()).items()}
        
    def thin(self,delta):
        """Take every delta non-unique samples."""
        c=ceil(cumsum([0]+self['weight'])/float(delta))
        ids=where(c[1:]>c[:-1])[0]
        weight=diff(c[[0]+list(ids+1)])
        t=self.sample(ids)
        t['weight']=weight
        return t
    
    def savecov(self,file,params=None):
        """Write the covariance to a file where the first line is specifies the parameter names."""
        if not params: params = self.params()
        with open(file,'w') as f:
            f.write("# "+" ".join(params)+"\n")
            savetxt(f,self.cov(params))
            
    def savechain(self,file,params=None):
        """Write the chain to a file where the first line is specifies the parameter names."""
        keys = ['lnl','weight']+list(params if params else self.params())
        with open(file,'w') as f:
            f.write("# "+" ".join(keys)+"\n")
            savetxt(f,self.matrix(keys))
            
    def plot(self,param):
        """Plot the value of a parameter as a function of sample number."""
        from matplotlib.pyplot import plot
        plot(cumsum(self['weight']),self[param])
        
    def like1d(self,p,**kwargs): 
        """Plots 1D likelihood contours for a parameter."""
        like1d(self[p],weights=self["weight"],**kwargs)
        
    def like2d(self,p1,p2,**kwargs): 
        """Plots 2D likelihood contours for a pair of parameters."""
        like2d(self[p1], self[p2], weights=self["weight"], **kwargs)
        
    def likegrid(self,**kwargs):
        """Plot several 2d likelihood contours."""
        if 'color' in kwargs: kwargs['colors']=[kwargs.pop('color')]
        likegrid(self,**kwargs)

        
        
class Chains(list):
    """A list of chains, probably from several MPI runs"""
    
    def burnin(self,n): 
        """Remove the first n samples from each chain."""
        return Chains(c.burnin(n) for c in self)
    
    def join(self): 
        """Combine the chains into one."""
        return Chain((k,hstack([c[k] for c in self])) for k in self[0].keys())
    
    def plot(self,param): 
        """Plot the value of a parameter as a function of sample number for each chain."""
        for c in self: c.plot(param)
    
def like2d(datx,daty,weights=None,
           nbins=15,which=[.68,.95],
           filled=True, color=None, cmap=None,
           ax=None, fig=None,
           **kwargs):
    
    from matplotlib.pyplot import figure, get_cmap
    from matplotlib.mlab import movavg
    from matplotlib.colors import LinearSegmentedColormap
    
    if ax is None: ax = (figure(0) if fig is None else (figure(fig) if isinstance(fig,int) else fig)).add_subplot(111)
    if weights is None: weights=ones(len(datx))
    if 'c' in kwargs: color=kwargs.pop('c')
    
    H,xe,ye = histogram2d(datx,daty,nbins,weights=weights)
    xem, yem = movavg(xe,2), movavg(ye,2)
    
    args = (xem,yem,transpose(H))
    kwargs = dict(levels=confint2d(H, sorted(which)[::-1]+[0]),**kwargs)
    if cmap is None: 
        cmap = {'b':'Blues',
                'g':'Greens',
                'r':'Reds',
                'orange':'Oranges',
                'grey':'Greys'}.get(color)
        if cmap is None: cmap = LinearSegmentedColormap.from_list(None,['w',color])
        else: cmap = get_cmap(cmap)
        
    if filled: ax.contourf(*args,cmap=cmap,**kwargs)
    ax.contour(*args,colors=color,**kwargs)
    
def like1d(dat,weights=None,
           nbins=30,range=None,maxed=True,
           ax=None, fig=None,
           **kw):
    from matplotlib.pyplot import figure
    from matplotlib.mlab import movavg
    if ax is None: ax = (figure(0) if fig is None else (figure(fig) if isinstance(fig,int) else fig)).add_subplot(111)
    if weights is None: weights=ones(len(dat))
    H, xe = histogram(dat,bins=nbins,weights=weights,normed=True,range=range)
    if maxed: H=H/max(H)
    xem=movavg(xe,2)
    ax.plot(xem,H,**kw)

def get_covariance(data,weights=None):
    if (weights==None): return cov(data.T)
    else:
        mean = sum(data.T*weights,axis=1)/sum(weights)
        zdata = data-mean
        return dot(zdata.T*weights,zdata)/(sum(weights)-1)


def likegrid(chains, params=None, 
             colors=['b','orange','k'], filled=True,
             nbins1d=30, nbins2d=20, 
             labels=None,
             fig=None,
             size=2,
             legend_loc=None):
    """
    Make a grid (aka "triangle plot") of 1- and 2-d likelihood contours. 
    
    Parameters
    ----------
    
    chains : 
        one or a list of `Chain` objects
    params, optional : 
        list of parameter names which to include of None for 
        all parameters (default: None)
    fig, optional :
        figure of figure number in which to plot (default: figure(0))
    colors, optional : 
        colors to cycle through for plotting
    filled, optional :
        whether to fill in the contours (default: True)
    labels, optional :
        list of names for a legend
    nbins1d, optional : 
        number of bins for 1d plots (default: 30)
    nbins2d, optional :
        number of bins for 2d plots (default: 20)
    """
    from matplotlib.pyplot import figure, Line2D
    fig = figure(0) if fig is None else (figure(fig) if isinstance(fig,int) else fig)
    if type(chains)!=list: chains=[chains]
    if params==None: params = sorted(reduce(lambda x,y: set(x)&set(y), [c.params() for c in chains]))
    if size is not None: fig.set_size_inches(*([size*len(params)]*2))
    colors=colors[:len(chains)]
    fig.subplots_adjust(hspace=0,wspace=0)

    c=chains[0]
    lims = [(max(min(c[p]),mean(c[p])-4*std(c[p])),min(max(c[p]),mean(c[p])+4*std(c[p]))) for p in params]
    ticks = [[t for t in ts if l[0]<=t<=l[1]] for (ts,l) in zip((c.mean(params)+c.std(params)*transpose([[-2,0,2]])).T ,lims)]
    n=len(params)
    for (i,p1) in enumerate(params):
        for (j,p2) in enumerate(params):
            if (i<=j):
                ax=fig.add_subplot(n,n,j*n+i+1)
                ax.set_xlim(*lims[i])
                ax.set_xticks(ticks[i])
                if (i==j): 
                    for (ch,col) in zip(chains,colors): 
                        if p1 in ch: ch.like1d(p1,nbins=nbins1d,color=col,ax=ax)
                    ax.set_yticks([])
                    
                elif (i<j): 
                    for (ch,col) in zip(chains,colors): 
                        if p1 in ch and p2 in ch: ch.like2d(p1,p2,filled=filled,nbins=nbins2d,color=col,ax=ax)
                    ax.set_ylim(*lims[j])
                    ax.set_yticks(ticks[j])
                        
                if i==0: 
                    ax.set_ylabel(p2)
                    ax.set_yticklabels(['%.3g'%t for t in ticks[j]])
                else: 
                    ax.set_yticklabels([])
                
                if j==n-1: 
                    ax.set_xlabel(p1)
                    ax.set_xticklabels(['%.3g'%t for t in ticks[i]])
                else: 
                    ax.set_xticklabels([])
                    
    fig.autofmt_xdate(rotation=90)
    
    if labels is not None:
        fig.legend([Line2D([0],[0],c=c) for c in colors],labels,fancybox=True,shadow=True,loc=legend_loc)


def confint2d(hist,which):
    """
    Return confidence levels in a histogram.
    
    Parameters
    ----------
    
    hist:
        A 2D histogram
    which: 
        A list of confidence levels, e.g. [.68, .95]
    """
    H=sort(hist.ravel())[::-1]
    sumH=sum(H)
    cdf=array([sum(H[H>x])/sumH for x in H])
    return interp(which,cdf,H)


def load_chain(path,paramnames=None):
    """
    Load a chain from a file or files in a variety of different formats.
    
    If path is a file, return a Chain object. The names of the parameters
    are expected either as a whitespace-separated comment on the first line of 
    the file (CosmoSlik style) or in a separate file 
    called <path>.paramnames (CosmoMC style). 
    
    If path is a directory, assumes it contains one file for each parameter (WMAP style)
    and gets the parameter name from the file name. 
    
    If path is a prefix such that there exists path_1, path_2, etc... 
    then returns a Chains object which is just a list of Chain objects 
    for each chain. 
    """
    
    if path.endswith('.ini'): 
        p = load_ini(path)
        return load_chain(p['output_file'])
    else:
        def load_one_chain(path):
            if os.path.isdir(path):
                if paramnames!=None: raise Exception("Can't specify custom parameter names if loading chain from a directory.")
                chain = {}
                for k in os.listdir(path):
                    try: chain[k]=loadtxt(os.path.join(path,k),usecols=[-1])
                    except: pass
                return Chain(chain)
            else:
                names = None
                if paramnames==None:
                    pnfiles = [os.path.join(os.path.dirname(path),f) for f in os.listdir(os.path.dirname(path)) if f.endswith('.paramnames') and os.path.basename(path).startswith(f[:-len('.paramnames')])]
                    if len(pnfiles)>1: raise Exception('Found multiple paramnames files for this chain; %s'%pnfiles)
                
                if paramnames or pnfiles:
                    with open(paramnames or pnfiles[0]) as f:
                        names = ['weight','lnl']+[line.split()[0] for line in f]
                        
                with open(path) as f:
                    if names==None: names = re.sub("#","",f.readline()).split()
                    try: data = loadtxt(f).T
                    except: data = [array([])]*len(names)
                    
                return Chain(zip(names,data))
        
        path = os.path.abspath(path)
        dir = os.path.dirname(path)
        files = [os.path.join(dir,f) for f in os.listdir('.' if dir=='' else dir) if re.match(os.path.basename(path)+'_[0-9]+',f) or f==os.path.basename(path)]
        if len(files)==1: return load_one_chain(files[0])
        elif len(files)>1: return Chains(filter(lambda c: c, (load_one_chain(f) for f in files)))
        else: raise IOError("File not found: "+path) 


