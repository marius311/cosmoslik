from numpy import array, loadtxt, dot, arange, diag, hstack, zeros
from scipy.linalg import cho_factor, cho_solve
from cosmoslik.plugins import Likelihood
from itertools import combinations_with_replacement
import os

class spt_r11(Likelihood):
    """
    
    ========================
    SPT Reichardt et al 2011
    ========================
    
    - Paper: `<http://adsabs.harvard.edu/abs/2011arXiv1111.0932R>`_
    - Data: `<http://pole.uchicago.edu/public/data/reichardt11/index.html>`_
    - CosmoSlik module by Marius Millea
    
    Usage
    =====
    
    To use this module add ``spt_r11`` to the list of ``likelihoods``
    
    
    Models
    ======
    
    This module required ``Models`` for 
    
    - cl_TT
    - Extra-galactic foregrounds
    
    
    Parameters
    ==========
    
    .. describe:: [spt_r11].a90
    .. describe:: [spt_r11].a150
    .. describe:: [spt_r11].a220
        
        These are the calibration factors at each frequency, defined so that 
        they multiply the theory spectrum. Calibration priors
        are included in the likelihood and are 1.75%, 1.6%, and 2.4% respectively.
        
        
    Plotting
    ========
    
    You can use the Inspector module to plot best-fit models and residuals. 
    ::
    
        import cosmoslik
        p = cosmoslik.inspect('param.ini')
        p['_likelihoods']['spt_r11'].plot(p=p,delta=delta)
    
    where ``delta`` is ``True`` if you want to plot residuals, otherwise ``False``.
    """



    def lnl(self, p, model):
            
        cl = self.get_cl_model(p,model)
            
        cl_vector = hstack([cl[spec_name] for spec_name in self.spec_names])
        
        dcl = self.spec_vector - cl_vector
        return dot(dcl,cho_solve(self.cov, dcl))/2 + self.lnl_calib(p)


    def lnl_calib(self,p):
        return sum(p.get(('spt_r11','a%s'%fr),1)-1**2/2/sig**2 \
                    for fr,sig in [('90',0175),('150',.016),('220',.024)])
        
        
    def get_cl_model(self,p,model):
        def get_cl(fr1,fr2): 
            calib = p.get(('spt_r11','a%s'%fr1),1)*p.get(('spt_r11','a%s'%fr2),1) 
            return calib * hstack([model['cl_TT'],zeros(10000)])[:self.lmax] + model['egfs']('cl_TT', lmax=self.lmax, freqs=(self.freq[fr1],self.freq[fr2]), fluxcut=self.fluxcut)
        def apply_windows(cl,windows): return array([dot(cl[self.windowrange],w[:,1]) for w in windows])
        return {spec_name:apply_windows(get_cl(*spec_name),windows) for (spec_name, windows) in self.windows.items()}
    
    
    def to_matrix(self,spec):
        return hstack([spec[spec_name] for spec_name in self.spec_names])
        
        
    def get_required_models(self, p):
        return ['cl_TT', 'egfs']
    
    
    def init(self, p):
        
        self.datadir     = os.path.join(os.path.dirname(__file__),'spt_multif_0809')
        
        self.spec_names  = [('90','90'),('90','150'),('90','220'),('150','150'),('150','220'),('220','220')]
        self.spec_vector = loadtxt(os.path.join(self.datadir,"spectrum_90_90x150_90x220_150_150x220_220.txt"))[:,1]
        self.spec        = dict(zip(self.spec_names,self.spec_vector.reshape(6,15)))
        self.cov         = cho_factor(loadtxt(os.path.join(self.datadir,"covariance_90_90x150_90x220_150_150x220_220.txt")).reshape(90,90))
        self.sigmas      = dict(zip(self.spec_names,diag(self.cov[0]).reshape(6,15)))
        self.windows     = dict(zip(self.spec_names,(lambda w: w.reshape(6,15,w.shape[1],2))(array([loadtxt(os.path.join(self.datadir,'window_%i'%i)) for i in range(1,91)]))))
        self.windowrange = (lambda w: slice(w[0,0],w[-1,0]+1))(self.windows.values()[0][0])
        self.ells        = {frs: array([dot(arange(self.windowrange.start,self.windowrange.stop),w[:,1]) for w in windows]) for frs,windows in self.windows.items()} 
        self.lmax        = self.windowrange.stop
        self.freq        = {'90' : {'dust':97.9,  'radio': 95.3,  'tsz':97.6},
                            '150': {'dust':153.8, 'radio': 150.2, 'tsz':152.9},
                            '220': {'dust':219.6, 'radio': 214.1, 'tsz':218.1}}
        self.fluxcut     = 6.4



    def plot(self, 
             fig=None, 
             cl=None, 
             p=None, 
             residuals=False, 
             show_data=True,
             show_model=True, 
             show_comps=False, 
             data_kw = {'c':'k'},
             model_kw = {'c':'k'},
             comps_kw = {},
             yscale='linear',
             ylim=None):
        
        from matplotlib.pyplot import figure
        if cl==None: cl = self.get_cl_model(p,p['_model'])
        if fig==None: fig=figure()
        fig.set_size_inches(15,15/1.6)
        fig.subplots_adjust(hspace=0,wspace=0)
        
        for ((i,fri),(j,frj)) in list(combinations_with_replacement(enumerate(['90','150','220']),2)):
            spec_name=(fri,frj)
            ax=fig.add_subplot(3,3,j*3+i+1)
            if residuals:
                if show_data: ax.errorbar(self.ells[spec_name],zeros(15),yerr=self.sigmas[spec_name],fmt='.',label='x'.join(spec_name),**data_kw)
                if show_model: ax.plot(self.ells[spec_name],cl[spec_name]-self.spec[spec_name],**model_kw)
                ax.set_ylim(-99,99)
            else:
                if show_data: ax.errorbar(self.ells[spec_name],self.spec[spec_name],yerr=self.sigmas[spec_name],fmt='.',label='x'.join(spec_name),**data_kw)
                if show_model: ax.plot(self.ells[spec_name],cl[spec_name],**model_kw)
                if show_comps: 
                    ax.plot(p['_model']['cl_TT'],c='b')
                    p['_model']['egfs']('cl_TT', lmax=self.lmax, freqs=(self.freq[fri],self.freq[frj]), fluxcut=self.fluxcut, plot=True, ax=ax, **comps_kw)
                ax.set_ylim(*(ylim or ((0,449) if yscale=='linear' else (1,1000))))
                ax.set_yscale(yscale)
            ax.set_xlim(1500,9500)
            if i==0: ax.set_ylabel(frj,size=16)
            else: ax.set_yticklabels([])
            if j==2: ax.set_xlabel(fri,size=16)
            else: ax.set_xticklabels([])

