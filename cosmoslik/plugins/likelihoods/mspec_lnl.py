import mspec as M
from numpy import dot, arange, diag, product
from scipy.linalg import cho_solve, cho_factor, cholesky
from cosmoslik.plugins import Likelihood
from itertools import combinations_with_replacement

class mspec_lnl(Likelihood):
    
    
    def get_required_models(self,p):
        return ['cl_TT']

    def init(self,p):
        if 'mspec' not in p: raise Exception('Expected an [mspec] section in the ini file.')
        
        self.mp = M.read_Mspec_ini(p['mspec'])
        
        self.signal = M.load_signal(self.mp)
        
        if 'rescale' in p['mspec']: self.signal = self.signal.rescaled(p['mspec'].get('rescale',1))
        if p['mspec'].get('to_dl'): self.signal = self.signal.dl()
        
        processed_signal = self.process_signal(p, self.signal, do_calib=False)

        self.lrange = self.mp['lrange']
        if isinstance(self.lrange,tuple): 
            self.lrange = {k:self.lrange for k in processed_signal.get_spectra()}
        
        self.signal_matrix_cov = processed_signal.get_as_matrix(lrange=self.lrange)[1]
        
        self.signal_matrix_cov = cholesky(self.signal_matrix_cov), False
        
        self.fluxcut = self.mp['fluxcut']
        self.eff_fr = self.mp['eff_fr']
        self.lmax = max([u for (_,u) in self.lrange.values()])
        
        
    def process_signal(self, p, sig=None, do_calib=True, keep_cov=True):
        """Get processed signal, appyling calibration, doing linear combination, etc.."""
        
        if sig is None: sig=self.signal
        
        sig = M.PowerSpectra(ells=sig.ells,
                             spectra=sig.spectra.copy(),
                             cov=(sig.cov.copy() if keep_cov else None),
                             binning=sig.binning)
        
        if do_calib:
            for (m1,m2) in sig.get_spectra():
                sig[m1,m2] = sig[m1,m2] * product([p['mspec'].get('calib',{}).get(m,1) for m in [m1,m2]])
        
        if 'cleaning' in p['mspec']: sig=sig.lincombo(p['mspec']['cleaning'])

        return sig

        
        
    def get_cl_model(self,p,model=None):
        """ 
        Build an Mspec PowerSpectra object which holds CMB + foreground C_ell's
        for all the required frequencies 
        """
        if model is None: model = p['_model']
        model_sig = M.PowerSpectra(ells=arange(self.lmax))
        if 'cleaning' in self.mp: 
            processed_spectra = M.utils.pairs({fr for coeffs in self.mp['cleaning'].values() for fr,w in coeffs if w!=0})
        else:
            processed_spectra = self.signal.get_spectra()
        for fr1,fr2 in processed_spectra:
            cl = model['cl_TT'][:self.lmax].copy()
            cl += model['egfs']('cl_TT',
                               fluxcut=min(self.fluxcut[fr1],self.fluxcut[fr2]),
                               freqs=(self.eff_fr[fr1],self.eff_fr[fr2]),
                               lmax=self.lmax)
            model_sig[(fr1,fr2)] = model_sig[(fr2,fr1)] = cl

        return self.process_signal(p,model_sig.binned(self.mp['binning']),do_calib=False)


    def plot(self,
             fig=None,
             cl=None, 
             p=None, 
             show_comps=False,
             show_model=True,
             yscale='log',
             ylim=None,
             residuals=False,
             data_color='k',
             model_color='k'):
        
        if cl==None: cl=self.get_cl_model(p, p['_model'])
        if fig==None: 
            from matplotlib.pyplot import figure
            fig=figure()
            
        processed_signal = self.process_signal(p)
        
        n=len(processed_signal.get_maps())

        fig.set_size_inches(6*n,6*n/1.6)
        fig.subplots_adjust(hspace=0,wspace=0)

        def slice_signal(sig,sl):      
            """Slice a signal according to an lrange"""
            return sig.sliced(sig.binning(slice(*sl)))
            
        for ((i,fri),(j,frj)) in combinations_with_replacement(enumerate(processed_signal.get_maps()),2):
            ax=fig.add_subplot(n,n,n*j+i+1)
            lrange = self.lrange[(fri,frj)]
            if residuals:
                slice_signal(processed_signal,lrange).diffed(slice_signal(cl,lrange)[fri,frj]).plot(ax=ax,which=[(fri,frj)],c=data_color)
            else:
                slice_signal(processed_signal,lrange).plot(ax=ax,which=[(fri,frj)],c=data_color)
                if show_model: 
                    cl.plot(ax=ax,which=[(fri,frj)],c=model_color)
                if show_comps:
                    ax.plot(p['_model']['cl_TT'],c='b')
                    p['_model']['egfs']('cl_TT',
                                        fluxcut=min(self.fluxcut[fri],self.fluxcut[frj]),
                                        freqs=(self.eff_fr[fri],self.eff_fr[frj]),
                                        lmax=self.lmax,
                                        plot=True,
                                        ax=ax)
                    
                ax.set_ylim(*(ylim or ((0,6999) if yscale=='linear' else (11,9999))))
                ax.set_yscale(yscale)
                ax.set_xlim(2,self.lmax-1)
                
                if n==1:
                    ax.set_title('%sx%s'%(fri,fri))
                else:
                    if i==0: ax.set_ylabel(frj,size=16)
                    else: ax.set_yticklabels([])
                    if j==n-1: ax.set_xlabel(fri,size=16)
                    else: ax.set_xticklabels([])



    
    def lnl(self,p,model):
        
        cl_model = self.get_cl_model(p, model)
        cl_model_matrix = cl_model.get_as_matrix(lrange=self.lrange).spec[:,1]
        signal_matrix_spec = self.process_signal(p, self.signal, keep_cov=False).get_as_matrix(lrange=self.lrange)[0][:,1]
        
        #Compute the likelihood  
        dcl = cl_model_matrix - signal_matrix_spec
        return dot(dcl,cho_solve(self.signal_matrix_cov,dcl))/2
