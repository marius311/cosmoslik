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
from numpy.linalg import inv
from itertools import takewhile, chain
from collections import defaultdict
import pickle
from functools import partial
from multiprocessing.pool import Pool
from numbers import Number
from functools import reduce
from subprocess import check_output, CalledProcessError, STDOUT


__all__ = ['Chain','Chains',
           'like1d','like2d','likegrid','likegrid1d','likepoints',
           'get_covariance', 'load_chain', 'load_cosmomc_chain']


class Chain(dict):
    """
    An MCMC chain. This is just a Python dictionary mapping parameter names
    to arrays of values, along with the special keys 'lnl' and 'weight'
    """
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
        for k,v in list(self.items()): self[k]=atleast_1d(v)
        if self and 'weight' not in self: self['weight']=ones(len(list(self.values())[0]))
        
    def copy(self):
        """Deep copy the chain so post-processing, etc... works right"""
        return Chain({k:v.copy() for k,v in self.items()})
        
    def params(self,non_numeric=False, fixed=False): 
        """
        Returns the parameters in this chain (i.e. the keys except 'lnl' and 'weight')
        Args:
            numeric: whether to include non-numeric parameters (default: False)
            fixed: whether to include parameters which don't vary (default: False)
        """
        return (set([k for k,v in list(self.items()) if 
                    (not k.startswith('_') 
                     and (non_numeric or (v.ndim==1 and issubclass(v.dtype.type,Number)))
                     and (fixed or not all(v==v[0]))
                     )])
                -set(["lnl","weight"]))
    
    def sample(self,s,keys=None): 
        """Return a sample or a range of samples depending on if s is an integer or a slice object."""
        return Chain((k,self[k][s]) for k in (keys if keys else list(self.keys())))

    def iterrows(self):
        """Iterate over the samples in this chain."""
        for i in range(self.length()):
            yield {k:v[i] for k,v in list(self.items())}

    def matrix(self,params=None):
        """Return this chain as an nsamp * nparams matrix."""
        if params is None: params=self.params()
        if is_iter(params) and not isinstance(params,str):
            return vstack([self[p] for p in (params if params else self.params())]).T
        else:
            return self[params]
    
    def cov(self,params=None): 
        """Returns the covariance of the parameters (or some subset of them) in this chain."""
        return get_covariance(self.matrix(params), self["weight"])
    
    def corr(self,params=None):
        """Returns the correlation matrix of the parameters (or some subset of them) in this chain."""
        return get_correlation(self.matrix(params), self["weight"])

    def mean(self,params=None): 
        """Returns the mean of the parameters (or some subset of them) in this chain."""
        return average(self.matrix(params),axis=0,weights=self["weight"])
    
    def std(self,params=None): 
        """Returns the std of the parameters (or some subset of them) in this chain."""
        return sqrt(average((self.matrix(params)-self.mean(params))**2,axis=0,weights=self["weight"]))
    
    def skew(self,params=None):
        """Return skewness of one or more parameters. """
        return average((self.matrix(params)-self.mean(params))**3,axis=0,weights=self["weight"])/self.std(params)**3

    def confbound(self,param,level=0.95,bound=None):
        """
        Compute an upper or lower confidence bound.
        
        Args:
            param (string): name of the parameter
            level (float): confidence level
            bound ('upper', 'lower', or None): whether to compute an upper or
                lower confidence bound. If set to None, will guess which based
                on the skewness of the distribution (will be lower bound for
                positively skewed distributions)
        """       
        if bound==None:
            bound = 'upper' if self.skew(param)>0 else 'lower'
        
        if bound=='lower':
            level = 1-level
            
        H, x = histogram(self[param],weights=self['weight'],bins=1000)
        xc = (x[1:]+x[:-1])/2
        b = interp(level,cumsum(H)/float(H.sum()),xc)
        return (self[param].min(),b) if bound=='upper' else (b,self[param].max())
    
    def acceptance(self): 
        """Returns the acceptance ratio."""
        return 1./mean(self["weight"])
    
    def length(self,unique=True):
        """Returns the number of unique samples. Set unique=False to get total samples."""
        return (len if unique else sum)(self['weight'])
    
    def burnin(self,n):
        """Remove the first n non-unique samples from the beginning of the chain."""
        return self.sample(slice(sum(1 for _ in takewhile(lambda x: x<n, cumsum(self['weight']))),None))

    def postprocd(self,func,nthreads=1,pool=None):
        """
        Post-process some values into this chain. 
                
        Args:
            func : a function which accepts all the keys in the chain
                and returns a dictionary of new keys to add. `func` must accept *all* 
                keys in the chain, if there are ones you don't need, capture them 
                with **_ in its call signature, e.g. to add in a parameter 'b'
                which is 'a' squared, use postprocd(lambda a,**_: {'b':a**2})
            nthreads : the number of threads to use
            pool : any worker pool which has a pool.map function. 
               default: multiprocessing.Pool(nthreads)
               
        Returns:
            A new chain with the new values post-processed in.
            Does not alter the original chain. If for some rows in the
            chain `func` did not return all the keys, these will be filled
            in with `nan`. 
            
        Note:
            This repeatedly calls `func` on rows in the chain, so its very inneficient 
            if you already have a vectorized version of your post-processing function. 
            `postprocd` is mostly useful for slow non-vectorized post-processing functions, 
            allowing convenient use of the `nthreads` option to this function. 
            
            For the default implementation of `pool`, `func` must be picklable, 
            meaning it must be a module-level function. 
        """

        if pool is not None: _pool = pool
        elif nthreads!=1: _pool = Pool(nthreads)
        else: _pool = None

        mp=map if _pool is None else _pool.map

        try:
            dat = list(mp(partial(_postprocd_helper,func),self.iterrows()))
        finally:
            if pool is None and _pool is not None: _pool.terminate()
        c=self.copy()
        allkeys = set(chain(*[list(d.keys()) for d in dat]))
        c.update({k:array([d.get(k,nan) for d in dat]) for k in allkeys})
        return c

    def reweighted(self,func,nthreads=1,pool=None):
        """
        Reweight this chain.

        Args:
            func : a function which accepts all keys in the chain, and returns 
                a new weight for the step. `func` must accept *all* 
                keys, if there are ones you don't need, capture them 
                with **_ in its call signature, e.g. to add unit gaussian prior 
                on parameter 'a' use reweighted(lambda a,**_: exp(-a**2/2)
            nthreads : the number of threads to use
            pool : any worker pool which has a pool.map function. 
               default: multiprocessing.Pool(nthreads)
               
        Returns:
            A new chain, without altering the original chain
            
        Note:
            This repeatedly calls `func` on rows in the chain, so its very inneficient 
            compared to a vectorized version of your post-processing function. `postprocd` is mostly 
            useful for slow post-processing functions, allowing you to conveniently 
            use the `nthreads` option to this function. 
            
            For the default implementation of `pool`, `func` must be picklable, 
            meaning it must be a module-level function. 
        """
        return self.postprocd(partial(_reweighted_helper,func),nthreads=nthreads,pool=pool)
    
    
    def add_gauss_prior(self, params, mean, covstd, nthreads=1, pool=None):
        """
        Post-process a gaussian prior into the chain.
        
        Args:
            params - a parameter name, or a list of parameters
            mean - the mean (should be a list if params was a list)
            covstd - if params was a list, this should be a 2-d array holding the covariance
                if params was a single parameter, this should be the standard devation
                
        Returns:
            A new chain, without altering the original chain
        """

        c=self.copy()   
        dx = self.matrix(params) - mean
        if is_iter(params) and not isinstance(params,str):
            dlnl=sum(dot(dx,inv(covstd))*dx,axis=1)/2
        else:
            dlnl=dx**2/2/covstd**2
        
        c['weight'] *= exp(-dlnl)
        c['lnl'] += dlnl
        return c

    def best_fit(self):
        """Get the best fit sample."""
        return {k:v[0] for k,v in list(self.sample(self['lnl'].argmin()).items())}
        
    def thin(self,delta):
        """Take every delta non-unique samples."""
        c=ceil(cumsum(append(0,self["weight"]))/float(delta))
        ids=where(c[1:]>c[:-1])[0]
        weight=diff(c[[0]+list(ids+1)])
        t=self.sample(ids)
        t['weight']=weight
        return t
        
    def thinto(self,num):
        """Thin so we end up with `num` total samples"""
        return self.thin(self.length(unique=False)//num)
    
    def savecov(self,file,params=None):
        """Write the covariance to a file where the first line is specifies the parameter names."""
        if not params: params = self.params()
        with open(file,'wb') as f:
            f.write(("# "+" ".join(params)+"\n").encode())
            savetxt(f,self.cov(params))
            
    def savechain(self,file,params=None):
        """Write the chain to a file where the first line is specifies the parameter names."""
        keys = ['lnl','weight']+list(params if params else self.params())
        with open(file,'w') as f:
            f.write("# "+" ".join(keys)+"\n")
            savetxt(f,self.matrix(keys))
            

    def plot(self,param=None,ncols=5,fig=None,size=4):
        """Plot the value of a parameter as a function of sample number."""
        from matplotlib.pyplot import figure
        if fig is None: fig=figure()
        params = [param] if param is not None else self.params()
        nrows = len(self.params())/ncols+1
        fig.set_size_inches(ncols*size,nrows*size/1.6)
        for i,param in enumerate(params,1):
            ax=fig.add_subplot(nrows,ncols,i)
            ax.plot(cumsum(self['weight']),self[param])
            ax.set_title(param)
        fig.tight_layout()

        
    def like1d(self,p,**kwargs): 
        """
        Plots 1D likelihood contours for a parameter.
        See :func:`~cosmoslik.chains.like1d`
        """
        like1d(self[p],weights=self["weight"],**kwargs)
        
    def like2d(self,p1,p2,**kwargs): 
        """
        Plots 2D likelihood contours for a pair of parameters.
        See :func:`~cosmoslik.chains.like2d`
        """
        like2d(self[p1], self[p2], weights=self["weight"], **kwargs)
        
    def likegrid(self,**kwargs):
        """
        Make a grid (aka "triangle plot") of 1- and 2-d likelihood contours. 
        See :func:`~cosmoslik.chains.likegrid`
        """
        if 'color' in kwargs: kwargs['colors']=[kwargs.pop('color')]
        likegrid(self,**kwargs)

    def likepoints(self,*args,**kwargs):
        """
        Plots samples from the chain as colored points.
        See :func:`~cosmoslik.chains.likepoints`
        """

        return likepoints(self,*args,**kwargs)

    def likegrid1d(self,**kwargs):
        """
        Make a grid of 1-d likelihood contours. 
        See :func:`~cosmoslik.chains.likegrid1d`
        """
        likegrid1d(self,**kwargs)

    def join(self):
        """Does nothing since already one chain."""
        return self
        
        
    def __repr__(self):
        return self.__str__()
        
    def __str__(self):
        """Print a summary of the chain. """
        return chain_stats(self)
        
        
def chain_stats(chain):
    lines = []
    lines.append(object.__repr__(chain))
    params = chain.params(fixed=True)
    maxlen = max(12,max(len(p) for p in params))+4
    if isinstance(chain,Chains):
        lines.append(('{:>%i}:  {:}'%maxlen).format('# of chains',len(chain)))
        chain = chain.join()
    lines.append(('{:>%i}:  {:}'%maxlen).format('# of steps',chain.length()))
    lines.append(('{:>%i}:  {:.2f}'%maxlen).format('total weight',chain['weight'].sum()))
    lines.append(('{:>%i}:  {:.3g}'%maxlen).format('acceptance',chain.acceptance()))
    for p in sorted(params):
        lines.append(('{:>%i}:  {:>10.4g} ± {:.4g}'%maxlen).format(p,chain.mean(p),chain.std(p)))
    return '\n'.join(lines)

    
        
class Chains(list):
    """A list of chains, e.g. from a run of several parallel chains"""
    
    def burnin(self,n): 
        """Remove the first n samples from each chain."""

        return Chains(c.burnin(n) for c in self)
    
    def join(self): 
        """Combine the chains into one."""
        return Chain((k,concatenate([c[k] for c in self])) for k in list(self[0].keys()))
    
    def plot(self,param=None,fig=None,**kwargs): 
        """Plot the value of a parameter as a function of sample number for each chain."""
        from matplotlib.pyplot import figure
        if fig is None: fig=figure()
        for c in self: 
            if c: c.plot(param,fig=fig,**kwargs)
            
    def params(self,**kwargs):
        return self[0].params(**kwargs)
            
    def __repr__(self):
        return self.__str__()
        
    def __str__(self):
        """Print a summary of the chain. """
        return chain_stats(self)


def _postprocd_helper(func,kwargs):
    return func(**kwargs)
    
def _reweighted_helper(func,weight,**kwargs):
    return {'weight': weight * func(weight=weight,**kwargs)}
        
def likepoints(chain,p1,p2,pcolor,
               npoints=1000,cmap=None,nsig=3,clim=None,marker='.',markersize=10,
               ax=None,zorder=-1,cbar=True,cax=None,**kwargs):
    """
    Plot p1 vs. p2 as points colored by the value of pcolor.
    
    Args:
        p1,p2,pcolor : parameter names
        npoints : first thins the chain so this number of points are plotted
        cmap : a colormap (default: jet)
        nsig : map the range of the color map to +/- nsig
        ax : axes to use for plotting (default: current axes)
        cbar : whether to draw a colorbar
        cax : axes to use for colorbar (default: steal from ax)
        marker, markersize, zorder, **kwargs : passed to the plot() command
    """
    from matplotlib.pyplot import get_cmap, cm, gca, sca, colorbar
    from matplotlib import colors, colorbar
    if cmap is None: cmap=get_cmap('jet')
    if ax is None: ax=gca()
    mu,sig = chain.mean(pcolor), chain.std(pcolor)
    for s in chain.thin(int(sum(chain['weight'])/float(npoints))).iterrows():
        if clim is None: c=cmap((s[pcolor]-mu)/(2*nsig*sig) + 0.5)
        else: c = cmap((s[pcolor]-clim[0])/(clim[1]-clim[0]))
        ax.plot(s[p1],s[p2],color=c,markeredgecolor=c,marker=marker,markersize=markersize,zorder=-1,**kwargs)
        
    if cax is None: cax = colorbar.make_axes(ax)[0]
    if clim is None: cb = colorbar.ColorbarBase(ax=cax,  norm=colors.Normalize(vmin=mu-nsig*sig, vmax=mu+nsig*sig))
    else: cb = colorbar.ColorbarBase(ax=cax,  norm=colors.Normalize(vmin=clim[0], vmax=clim[1]))
    
    sca(ax)
    return ax,cax


def like2d(datx,daty,weights=None,
           nbins=15,which=[.68,.95],
           filled=True, color=None, cmap=None,
           ax=None,
           **kwargs):
    
    from matplotlib.pyplot import gca, get_cmap
    from matplotlib.mlab import movavg
    from matplotlib.colors import LinearSegmentedColormap
    
    if ax is None: ax = gca()
    if weights is None: weights=ones(len(datx))
    if color is None: color = kwargs.pop('c') if 'c' in kwargs else 'b' 
    
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
           nbins=30,ranges=None,maxed=True,
           ax=None,smooth=False,
           kde=True,
           zero_endpoints=False,
           filled=False,
           **kw):
    
    from matplotlib.pyplot import gca
    if ax is None: ax = gca()
    
    if kde:
        try:
            from getdist import MCSamples
        except ImportError as e:
            raise Exception("Plotting with kde, kde1d, or kde2d set to True requires package 'getdist'. Install this package or set to False.") from e
        
        d = MCSamples(samples=dat, weights=weights, names=['x'], ranges={'x':ranges or (None,None)}).get1DDensity(0)
        d.normalize('max' if maxed else 'integral')
        xem, H = d.x, d.P * (maxed or 1)
        
    else:
    
        from matplotlib.mlab import movavg
        H, xe = histogram(dat,bins=nbins,weights=weights,normed=True,range=ranges)
        xem=movavg(xe,2)
        
        if smooth:
            from scipy.interpolate import PchipInterpolator
            itp = PchipInterpolator(xem,H)
            xem = linspace(xem.min(),xem.max(),100)
            H = itp(xem) 
        
        if maxed: H = H/max(H) * (maxed or 1)
        
    if zero_endpoints:
        xem = hstack([[xem[0]],xem,[xem[-1]]])
        H = hstack([[0],H,[0]])
    
    if filled:
        ax.fill_between(xem,H,**kw)
    else:
        ax.plot(xem,H,**kw)

def get_correlation(data,weights=None):
    cv = get_covariance(data,weights)
    n,n = cv.shape
    for i in range(n): 
        std=sqrt(cv[i,i])
        cv[i,:]/=std
        cv[:,i]/=std
    return cv

def get_covariance(data,weights=None):
    if (weights is None): return cov(data.T)
    else:
        mean = sum(data.T*weights,axis=1)/sum(weights)
        zdata = data-mean
        return dot(zdata.T*weights,zdata)/(sum(weights)-1)


def likegrid(chains, params=None, 
             lims=None, ticks=None, nticks=4,
             nsig=4,
             spacing=0.05,
             xtick_rotation=30,
             colors=None, filled=True,
             nbins1d=30, smooth1d=False, kde1d=True,
             nbins2d=20,
             labels=None,
             fig=None,
             size=2,
             legend_loc=None,
             param_name_mapping=None,
             param_label_size=None):
    """
    Make a grid (aka "triangle plot") of 1- and 2-d likelihood contours. 
    
    Parameters
    ----------
    
    chains : 
        one or a list of `Chain` objects
        
    default_chain, optional :
        the chain used to get default parameters names, axes limits, and ticks 
        either an index into chains or a `Chain` object (default: chains[0])
        
    params, optional : 
        list of parameter names which to show 
        (default: all parameters from default_chain)
        
    lims, optional :
        a dictionary mapping parameter names to (min,max) axes limits
        (default: +/- 4 sigma from default_chain)
        
    ticks, optional :
        a dictionary mapping parameter names to list of [ticks]
        (default: automatically picks `nticks`)

    nticks, optional :
        roughly how many ticks per axes (default: 5)
        
    xtick_rotation, optional :
        numbers of degrees to rotate the xticks by (default: 30)

    spacing, optional :
        space in between plots as a fraction of figure width (default: 0.05)

    fig, optional :
        figure of figure number in which to plot (default: figure(0))
        
    size, optional :
        size in inches of one plot (default: 2)
        
    colors, optional : 
        colors to cycle through for plotting
        
    filled, optional :
        whether to fill in the contours (default: True)
        
    labels, optional :
        list of names for a legend
        
    legend_loc, optional :
        (x,y) location of the legend (coordinates scaled to [0,1]) 
        
    nbins1d, optional : 
        number (or len(chains) length list) of bins for 1d plots (default: 30)
        
    nbins2d, optional :
        number (or len(chains) length list) of bins for 2d plots (default: 20)
    """
    from matplotlib.pyplot import figure, Line2D, xticks
    from matplotlib.ticker import MaxNLocator
    fig = figure(0) if fig is None else (figure(fig) if isinstance(fig,int) else fig)
    if type(chains)!=list: chains=[chains]
    if params==None: params = sorted(reduce(lambda x,y: set(x)&set(y), [c.params() for c in chains]))
    if param_name_mapping is None: param_name_mapping = {}
    if size is not None: fig.set_size_inches(*([size*len(params)]*2))
    if colors is None: colors=['b','orange','k','m','cyan']
    if not isinstance(nbins2d,list): nbins2d = [nbins2d]*len(chains)
    if not isinstance(nbins1d,list): nbins1d = [nbins1d]*len(chains)
    fig.subplots_adjust(hspace=spacing,wspace=spacing)
    
    if lims is None: lims = {}
    lims = {p:(lims[p] if p in lims 
               else (min(max(min(c[p]),c.mean(p)-nsig*c.std(p)) for c in chains if p in c.params()),
                     max(min(max(c[p]),c.mean(p)+nsig*c.std(p)) for c in chains if p in c.params()))) 
            for p in params}
    if ticks is None: ticks = {}
    if isinstance(nticks,int): nticks={p:nticks for p in params}

    n=len(params)
    for (i,p1) in enumerate(params):
        for (j,p2) in enumerate(params):
            if (i<=j):
                ax=fig.add_subplot(n,n,j*n+i+1)
                ax.xaxis.set_major_locator(MaxNLocator(nticks.get(p1,5)))
                ax.yaxis.set_major_locator(MaxNLocator(nticks.get(p2,5)))
                ax.set_xlim(*lims[p1])
                if (i==j): 
                    for (ch,col,nbins) in zip(chains,colors,nbins1d): 
                        if p1 in ch: ch.like1d(p1,nbins=nbins,color=col,ax=ax,smooth=smooth1d,kde=kde1d)
                    ax.set_yticklabels([])
                    
                elif (i<j): 
                    for (ch,col,nbins) in zip(chains,colors,nbins2d): 
                        if p1 in ch and p2 in ch: ch.like2d(p1,p2,filled=filled,nbins=nbins,color=col,ax=ax)
                    if p2 in ticks: ax.set_yticks(ticks[p2])
                    ax.set_ylim(*lims[p2])
                        
                if i==0: 
                    ax.set_ylabel(param_name_mapping.get(p2,p2),size=param_label_size)
                else: 
                    ax.set_yticklabels([])
                
                if j==n-1: 
                    ax.set_xlabel(param_name_mapping.get(p1,p1),size=param_label_size)
                    xticks(rotation=xtick_rotation)
                else: 
                    ax.set_xticklabels([])
                    
    if labels is not None:
        fig.legend([Line2D([0],[0],c=c,lw=2) for c in colors],labels,
                   fancybox=True,shadow=False,
                   loc='upper right', bbox_to_anchor=(legend_loc or (0.8,0.8)))


from collections import Iterable
import operator as op

def likegrid1d(chains, 
               params='all',
               lims=None, 
               ticks=None,
               nticks=4,
               nsig=3,
               colors=None,
               nbins1d=30,
               smooth1d=False,
               kde1d=True,
               labels=None,
               fig=None,
               size=2,
               aspect=1,
               legend_loc=None,
               linewidth=1,
               param_name_mapping=None,
               param_label_size=None,
               tick_label_size=None,
               titley=1,
               ncol=4,
               axes=None):
    """
    Make a grid of 1-d likelihood contours.
   
    Arguments:
    ----------
   
    chains :
        one or a list of `Chain` objects
       
    default_chain, optional :
        the chain used to get default parameters names, axes limits, and ticks
        either an index into chains or a `Chain` object (default: chains[0])
       
    params, optional :
        list of parameter names which to show
        can also be 'all' or 'common' which does the union/intersection of
        the params in all the chains
       
    lims, optional :
        a dictionary mapping parameter names to (min,max) axes limits
        (default: +/- 4 sigma from default_chain)
       
    ticks, optional :
        a dictionary giving a list of ticks for each parameter
        
    nticks, optional :
        roughly how many x ticks to show. can be dictionary to 
        specify each parameter separately. (default: 4)
       
    fig, optional :
        figure of figure number in which to plot (default: new figure)

    ncol, optional :
        the number of colunms (default: 4)

    axes, optional :
        an array of axes into which to plot. if this is provided, fig and ncol 
        are ignored. must have len(axes) >= len(params). 
       
    size, optional :
        size in inches of one plot (default: 2)

    aspect, optional :
        aspect ratio (default: 1)

       
    colors, optional :
        colors to cycle through for plotting
       
    filled, optional :
        whether to fill in the contours (default: True)
       
    labels, optional :
        list of names for a legend
       
    legend_loc, optional :
        (x,y) location of the legend (coordinates scaled to [0,1])
       
    nbins1d, optional :
        number of bins for 1d plots (default: 30)
       
    nbins2d, optional :
        number of bins for 2d plots (default: 20)


    """
    from matplotlib.pyplot import gcf, Line2D
    from matplotlib.ticker import AutoMinorLocator, ScalarFormatter, MaxNLocator

    if type(chains)!=list: chains=[chains]
        
    if params in ['all','common']: 
        params = sorted(reduce(lambda x,y: (op.__or__ if params=='all' else op.__and__)(set(x),set(y)), [c.params() for c in chains]))
    elif not isinstance(params,Iterable):
        raise ValueError("params should be iterable or 'all' or 'common'")
                         
    if param_name_mapping is None: param_name_mapping = {}
    nrow = int(floor(len(params)/ncol))+1
    if fig is None:
        fig = gcf()
        fig.subplots_adjust(hspace=0.4,wspace=0.1)#,bottom=0, top=1, left=0, right=1)
        if size is not None: 
            fig.set_size_inches(size*ncol,size*nrow/aspect)
                            
    if colors is None: 
        colors = ['b','orange','k','m','cyan']
       
    if lims is None: lims = {}
    lims = {p:(lims[p] if p in lims 
               else (min(max(min(c[p]),c.mean(p)-nsig*c.std(p)) for c in chains if p in c.params()),
                     max(min(max(c[p]),c.mean(p)+nsig*c.std(p)) for c in chains if p in c.params()))) 
            for p in params}
    
    for (i,p1) in enumerate(params,1):
        ax=fig.add_subplot(nrow,ncol,i)
        if ticks is not None and p1 in ticks:
            ax.set_xticks(ticks[p1])
        for (ch,col) in zip(chains,colors):
            if p1 in ch: ch.like1d(p1,nbins=nbins1d,color=col,ax=ax,linewidth=linewidth,smooth=smooth1d,kde=kde1d)
        ax.set_yticks([])
        ax.set_xlim(lims[p1])
        ax.set_ylim(0,1)
        ax.set_title(param_name_mapping.get(p1,p1),size=param_label_size,y=titley)
        ax.tick_params(labelsize=tick_label_size)
        if ticks and p1 in ticks:
            ax.set_xticks(ticks[p1])
        else:
            ax.xaxis.set_major_locator(MaxNLocator(nbins=nticks.get(p1,4) if isinstance(nticks,dict) else nticks))
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.xaxis.set_major_formatter(ScalarFormatter(useOffset=False))

    if labels is not None:
        lg = fig.legend([Line2D([0],[0],c=c,linewidth=2) for c in colors],labels,fancybox=True,shadow=False,
                   loc='upper center', bbox_to_anchor=(legend_loc or (0.5,1.1)))
        return lg
        
        
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




def load_cosmomc_chain(path,paramnames=None):
    """
    Load a chain from a file or files in a variety of different formats.
    
    If ``path`` is a CosmoSlik ini, the read the ``output_file`` key from the 
    ini and load that chain.
    
    If ``path`` is a file, return a :class:`~cosmoslik.chains.Chain` object. The names of the parameters
    are expected either as a whitespace-separated comment on the first line of 
    the file (CosmoSlik style) or in a separate file 
    called <path>.paramnames (CosmoMC style). 
    
    If ``path`` is a directory, assumes it contains one file for each parameter (WMAP style)
    and gets the parameter name from the file name. 
    
    If ``path`` is a prefix such that there exists <path>_1, <path>_2, etc... 
    then returns a :class:`~cosmoslik.chains.Chains` object which is a list of chains.
    
    """
    
    if path.endswith('.ini'): 
        p = load_ini(path)
        return load_chain(p['output_file'])
    else:
        def load_one_chain(path):
            
            nonlocal paramnames
            
            if os.path.isdir(path):
                
                if paramnames is not None: raise Exception("Can't specify custom parameter names if loading chain from a directory.")
                chain = {}
                for k in os.listdir(path):
                    try: chain[k]=loadtxt(os.path.join(path,k),usecols=[-1])
                    except: pass
                return Chain(chain)
                
            else:
                
                # try automatically finding the corresponding *.paramnames file
                if paramnames is None:
                    pnfiles = [os.path.join(os.path.dirname(path),f) for f in os.listdir(os.path.dirname(path)) if f.endswith('.paramnames') and os.path.basename(path).startswith(f[:-len('.paramnames')])]
                    if len(pnfiles)>1: 
                        raise Exception('Found multiple paramnames files for this chain; %s'%pnfiles)
                    elif len(pnfiles)==1: 
                        paramnames = pnfiles[0]
                
                # if we have a *.paramnames file at this point, load it
                if isinstance(paramnames,str):
                    with open(paramnames) as f:
                        paramnames = ['weight','lnl']+[line.split()[0] for line in f]
                        
                with open(path) as f:
                    
                    # if still no *.paramnames, look for names inside a comment on the first line
                    if paramnames is None: 
                        line = f.readline()
                        if not line.startswith("#"):
                            raise Exception("Couldn't find any paramnames. Specify paramnames=... by hand.")
                        paramnames = re.sub("#","",line).split()
                            
                    try:
                        data = loadtxt(f).T
                    except: 
                        return None
                    else:
                        return Chain(list(zip(paramnames,data)))
                
                    
        
        path = os.path.abspath(path)
        dir = os.path.dirname(path)
        files = [os.path.join(dir,f) for f in os.listdir('.' if dir=='' else dir) if re.match(os.path.basename(path)+'_[0-9]+',f) or f==os.path.basename(path)]
        if len(files)==1: return load_one_chain(files[0])
        elif len(files)>1: return Chains([c for c in (load_one_chain(f) for f in files) if c])
        else: raise IOError("File not found: "+path) 


def is_iter(x):
    try: 
        iter(x)
        return True
    except: 
        return False


def combine_covs(*covs):
    """
    Combines a bunch of covariances into a single array. If a parameter appears
    in multiple covariances, the covariance appearing last in the list is used. 
    
    Each cov can be:
        * tuple of ([names...], 2-d array)
        * filename (file's first line is "#paramname1 paramname2 ..." and next lines are array)
        * `Chain` object (will call its cov() function)
        * {k:std(v)...}
    
    Returns:
        Tuple of ([names...], 2-d array)
    """
    def to_array(cov):
        if isinstance(cov,tuple): 
            return cov
        elif isinstance(cov,str):
            with open(cov) as f:
                return re.sub("#","",f.readline()).split(), loadtxt(f)
        elif isinstance(cov,Chain):
            return cov.params(), cov.cov()
        elif isinstance(cov,dict):
            return [k for k in cov], diag([v**2 for v in list(cov.values())])
        else:
            raise ValueError("Unrecognized covariance data type.")
        
    covs = [to_array(cov) for cov in covs]
    allnames = list(chain(*[n for n,_ in covs]))
    
    allcov = zeros((len(allnames),len(allnames)))
    for (names,cov) in covs:
        idxs = [allnames.index(n) for n in names]
        for i in idxs: allcov[i,:] = allcov[:,i] = 0
        allcov[ix_(idxs,idxs)] = cov
    
    return allnames, allcov


def load_chain(filename, repack=False):
    """
    Load a chain produced by a compatible CosmoSlik sampler like
    metropolis_hastings or emcee.
    
    
    Parameters:
    -----------
    repack : 
        If the chain file is not currently open (i.e. the chain is not currently
        running), and if the chain is stored in chunks as output from an MPI run,
        then overwrite the file with a more efficiently stored version which
        will be faster to load the next time. 
    
    """
    with open(filename, 'rb') as f:
        c = pickle.load(f)
        if isinstance(c,(Chain,Chains)): 
            return c
        else:
            names = [n.decode() if isinstance(n,bytes) else n for n in c]
            dat = []
            while True:
                try:
                    dat.append(pickle.load(f,encoding="latin1"))
                except:
                    break
            ii=set(i for i,_ in dat)

            if dat[0][1].dtype.kind=='V':
                c = Chains([Chain({n:concatenate([d['f%i'%k] for j,d in dat if i==j]) for k,n in enumerate(names)}) for i in ii])
            else:
                c = Chains([Chain(dict(list(zip(names,vstack([d for j,d in dat if i==j]).T)))) for i in ii])
                
    if repack:
        try:
            open_files = check_output(["lsof","--",filename], stderr=STDOUT).splitlines()[1:]
        except CalledProcessError as e:
            if e.returncode == 1 and e.output==b'':
                open_files = []
            else:
                return c
            
        if not any([l.split()[3].decode()[-1] in "uw" for l in open_files]):
            with open(filename,'wb') as f:
                pickle.dump(c,f)
                print("Repacked: "+filename)
                
    return c
