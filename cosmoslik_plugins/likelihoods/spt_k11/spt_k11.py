from numpy import array, fromstring, loadtxt, dot, arange, diag
from scipy.linalg import cho_factor, cho_solve
from cosmoslik.plugins import Likelihood
from itertools import takewhile
import os

class spt_k11(Likelihood):

    def lnl(self, p, model):
            
        cl = self.get_cl_model(model)
            
        if p.get('diagnostic',False):
            from matplotlib.pyplot import ion, figure, draw, cla
            ion()
            cla()
            ax = figure().add_subplot(111)
            self.plot(ax,'cl_TT',cl)
            ax.set_yscale('log')
            ax.set_ylim(10,6e3)
            draw()
            
        #Apply windows and calculate likelihood
        dcl = self.spec-cl
        return dot(dcl,cho_solve(self.sigma, dcl))/2


    def plot(self,p, cl=None,
             ax=None,fig=None,
             residuals=False,
             show_comps=True,
             comps_kw={}):
        
        if cl==None: cl = self.get_cl_model(p['_model'])
        
        if ax is None:
            if fig==None: 
                from matplotlib.pyplot import figure
                fig=figure()
            ax=fig.add_subplot(111)
            
        if not residuals:
            ax.errorbar(self.ells,self.spec,yerr=diag(self.sigma[0]),fmt='.',label='SPT K11')
            ax.plot(self.ells,cl)
            if show_comps: 
                ax.plot(p['_model']['cl_TT'],c='b')
                p['_model']['egfs']('cl_TT', lmax=self.lmax, freqs=(self.freq,self.freq), fluxcut=self.fluxcut, plot=True, ax=ax, **comps_kw)

        else:
            ax.errorbar(self.ells,self.spec-cl,yerr=diag(self.sigma[0]),fmt='.',label='SPT K11')
            ax.plot([self.ells[0],self.ells[-1]],[0]*2)


    def get_cl_model(self,model):
        #Get CMB + foreground model
        cl = model['cl_TT'][:self.lmax] + \
             model['egfs']('cl_TT', lmax=self.lmax, freqs=(self.freq,self.freq), fluxcut=self.fluxcut)
        
        #Apply window functions
        return array([dot(cl[self.windowrange],w) for w in self.windows])
    
    def get_required_models(self, p):
        return ['cl_TT', 'egfs']

    def init(self, p):
        
        self.datadir = os.path.join(os.path.dirname(__file__),'bandpowers')
        
        #Load spectrum and covariance
        with open(os.path.join(self.datadir,'Spectrum_spt20082009.newdat')) as f:
            while 'TT' not in f.readline(): pass
            self.spec=array([fromstring(f.readline(),sep=' ')[1] for _ in range(47)])
            self.sigma=array([fromstring(f.readline(),sep=' ') for _ in range(94)])[47:]
            
        #Load windows
        self.windows = [loadtxt(os.path.join(self.datadir,'windows','window_0809','window_%i'%i))[:,1] for i in range(1,48)]
        
        if ('spt_k11','lmin') in p:
            bmin = sum(1 for _ in takewhile(lambda x: x<p['spt_k11','lmin'], (sum(1 for _ in takewhile(lambda x: abs(x)<.001,w) ) for w in self.windows)))
            self.spec = self.spec[bmin:]
            self.sigma = self.sigma[bmin:,bmin:]
            self.windows = self.windows[bmin:]
        
        self.sigma = cho_factor(self.sigma)
        
        self.windowrange = (lambda x: slice(min(x),max(x)+1))(loadtxt(os.path.join(self.datadir,'windows','window_0809','window_1'))[:,0])
        self.lmax = self.windowrange.stop
        self.ells = array([dot(arange(10000)[self.windowrange],w) for w in self.windows])

        self.freq = {'dust':154, 'radio': 151, 'tsz':153}
        self.fluxcut = 50
        
    




