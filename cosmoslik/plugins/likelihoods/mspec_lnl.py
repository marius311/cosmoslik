import mspec as M
from mspec.utils import pairs
from numpy import dot, arange, product, zeros
from scipy.linalg import cho_solve, cholesky
from cosmoslik.plugins import Likelihood

class mspec_lnl(Likelihood):
    
    
    def get_required_models(self,p):
        return ['cl_TT']

    def init(self,p):
        if 'mspec' not in p: raise Exception('Expected an [mspec] section in the ini file.')
        
        self.mp = M.read_Mspec_ini(p['mspec'])
        
        self.signal = M.load_signal(self.mp)
        
        if 'rescale' in p['mspec']: self.signal = self.signal.rescaled(p['mspec'].get('rescale',1))
        if p['mspec'].get('to_dl'): self.signal = self.signal.dl()
        
        self.cleaning = self.mp['cleaning'] if 'cleaning' in self.mp else {m:[(m,1)] for m in self.signal.get_maps()}
        
        self.per_freq_egfs = M.SymmetricTensorDict(self.mp.get('per_freq_egfs',{}))
        
        processed_signal = self.process_signal(p, self.signal, do_calib=False)

        self.lrange = self.mp['lrange']
        if isinstance(self.lrange,tuple): 
            self.lrange = {k:self.lrange for k in processed_signal.get_spectra()}
        
        self.signal_matrix_cov = processed_signal.get_as_matrix(lrange=self.lrange)[1]
        
        self.signal_matrix_cov = cholesky(self.signal_matrix_cov), False
        
        self.fluxcut = self.mp.get('fluxcut')
        self.eff_fr = self.mp.get('eff_fr')
        self.lmax = max([u for (_,u) in self.lrange.values()])
        
        
    def process_signal(self, p, sig=None, do_calib=True, keep_cov=True):
        """Get processed signal, appyling calibration, doing linear combination, etc.."""
        
        if sig is None: sig=self.signal
        
        sig = M.PowerSpectra(ells=sig.ells,
                             spectra=sig.spectra.copy(),
                             cov=(sig.cov.copy() if keep_cov else None),
                             binning=sig.binning)
        
        def calib(sig):
            if do_calib:
                for (m1,m2) in sig.get_spectra():
                    sig[m1,m2] = sig[m1,m2] * product([p['mspec'].get('calib',{}).get(m,1) for m in [m1,m2]])
        
         
        calib(sig) #apply calibration to frequency PS
        sig=sig.lincombo(self.cleaning)
        calib(sig) #apply calibration to cleaned PS
            
        return sig

    def get_cl_model(self,p,
                     model=None,
                     get_cmb=True,
                     get_egfs=True,
                     egfs_kwargs=None):
        """ 
        Build an Mspec PowerSpectra object which holds CMB + foreground C_ell's
        for all the required frequencies 
        """
        if model is None: model = p['_model']
        model_sig = M.PowerSpectra(ells=arange(self.lmax))
            
        in_spectra = pairs({m for coeffs in self.cleaning.values() for m,_ in coeffs})
        out_spectra = pairs(self.cleaning.keys())
            
        def add_in_components(model_sig,spectra, default_egfs=True, get_cmb=True):
            for fr1,fr2 in spectra:
                cl = model_sig.spectra.setdefault((fr1,fr2),zeros(self.lmax))
                if get_cmb: 
                    cl += model['cl_TT'][:self.lmax].copy()
                if get_egfs and (default_egfs or (fr1,fr2) in self.per_freq_egfs):
                    cl += model['egfs']('cl_TT',
                                        p_egfs = p[self.per_freq_egfs[(fr1,fr2)]] if (fr1,fr2) in self.per_freq_egfs else None,
                                        fluxcut=min(self.fluxcut[fr1],self.fluxcut[fr2]) if self.fluxcut else None,
                                        freqs=(self.eff_fr[fr1],self.eff_fr[fr2]) if self.eff_fr else None,
                                        lmax=self.lmax)
                    cl += p.get(('mspec','gal','amp_%s_%s'%tuple(sorted([fr1,fr2]))),0) * (arange(self.lmax)/float(p.get(('mspec','gal','norm_ell'),3000)))**p.get(('mspec','gal','tilt'),0)

                    
            return model_sig

        if self.per_freq_egfs:
            model_sig = add_in_components(model_sig, in_spectra, get_cmb=get_cmb, default_egfs=False)
            model_sig = self.process_signal(p,model_sig,do_calib=False)
            model_sig = add_in_components(model_sig, out_spectra, get_cmb=False, default_egfs=False)
            return model_sig.binned(self.mp['binning'])
        else:
            return self.process_signal(p,add_in_components(model_sig, in_spectra, get_cmb=get_cmb).binned(self.mp['binning']),do_calib=False)


    def plot(self,
             fig=None,
             p=None, 
             show_comps=False,
             show_model=True,
             yscale='log',
             ylim=None,
             residuals=False,
             data_color='k',
             model_color='k'):
        
        model_total = self.get_cl_model(p)
        if show_comps: 
            model_cmb = self.get_cl_model(p, get_cmb=True, get_egfs=False)
            model_egfs = self.get_cl_model(p, get_cmb=False, get_egfs=True)
        
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
            
        for ((i,fri),(j,frj)) in pairs(enumerate(processed_signal.get_maps())):
            ax=fig.add_subplot(n,n,n*j+i+1)
            lrange = self.lrange[(fri,frj)]
            if residuals:
                slice_signal(processed_signal,lrange).diffed(slice_signal(model_total,lrange)[fri,frj]).plot(ax=ax,which=[(fri,frj)],c=data_color)
                ax.set_ylim(*(ylim or (-39,39)))
            else:
                slice_signal(processed_signal,lrange).plot(ax=ax,which=[(fri,frj)],c=data_color)
                if show_model: 
                    slice_signal(model_total,lrange).plot(ax=ax,which=[(fri,frj)],c=model_color)
                if show_comps:
                    slice_signal(model_cmb,lrange).plot(ax=ax,which=[(fri,frj)],c='b')
                    slice_signal(model_egfs,lrange).plot(ax=ax,which=[(fri,frj)],c='orange')
                    
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
