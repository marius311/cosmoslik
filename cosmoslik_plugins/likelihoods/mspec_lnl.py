import mspec as M
from mspec.utils import pairs
from numpy import dot, arange, product, zeros, load, exp, vstack, hstack, ones, loadtxt, array
from scipy.linalg import cho_solve, cholesky
from cosmoslik.plugins import Likelihood

class mspec_lnl(Likelihood):
    
    def get_required_models(self,p):
        return ['cl_TT']

    def init(self,p):
        if 'mspec' not in p: raise Exception('Expected an [mspec] section in the ini file.')
        
        if 'mspec_ini' in p['mspec']: 
            self.signal = M.load_signal(p['mspec','mspec_ini'])
            self.mp = p['mspec']
            self.mp.update(M.read_Mspec_ini(p['mspec','mspec_ini']))
        else: 
            self.mp = M.read_Mspec_ini(p['mspec'])
            self.signal = M.load_signal(self.mp)
        
        if 'rescale' in self.mp: self.signal = self.signal.rescaled(self.mp['rescale'])
        if self.mp.get('to_dl'): self.signal = self.signal.dl()
        
        self.cleaning = self.mp['cleaning'] if 'cleaning' in self.mp else {m:[(m,1)] for m in self.signal.get_maps()}
        
        self.beamcov = None
        if 'beam' in self.mp:
            self.beampca = M.SymmetricTensorDict()
            for k,v in self.mp['beam'].items():
                self.beampca[tuple(k.split('_'))] = load(v['file'])
            if 'beamcov' in self.mp:
                with open(self.mp['beamcov']) as f:
                    self.beamcov = (f.readline().replace('#','').split(),(cholesky(loadtxt(f)),False))
                
        self.per_freq_egfs = M.SymmetricTensorDict(self.mp.get('per_freq_egfs',{}))
        
        self.fluxcut = self.mp.get('fluxcut')
        self.eff_fr = self.mp.get('eff_fr')
        self.lrange = self.mp['lrange']
        self.lmax = self.lrange[1] if isinstance(self.lrange,tuple) else max([u for (_,u) in self.lrange.values()])

        processed_signal = self.process_signal(p, self.signal, do_calib=False)

        if isinstance(self.lrange,tuple): 
            self.lrange = {k:self.lrange for k in processed_signal.get_spectra()}
        
        self.signal_matrix_cov = processed_signal.get_as_matrix(lrange=self.lrange)[1]
        self.signal_matrix_cov = cholesky(self.signal_matrix_cov), False
        
        
        
    def process_signal(self, p, sig=None, do_calib=True, do_beam=True, keep_cov=True):
        """Get processed signal, appyling calibration, doing linear combination, etc.."""
        
        if sig is None: sig=self.signal
        
        sig = sig.deepcopy()
        
        def calib(sig):
            if do_calib:
                for (m1,m2) in sig.get_spectra():
                    sig[m1,m2] = sig[m1,m2] * product([p['mspec'].get('calib',{}).get(m,1) for m in [m1,m2]])
        
         
        calib(sig) #apply calibration to frequency PS
        
        if do_beam and 'beam' in self.mp :
            for ps,dbl in self.get_beam_correction(p).items():
                dbl = sig.binning(dbl)
                nl = min(sig[ps].size,dbl.size)
                sig[ps][:nl] = sig[ps][:nl] * dbl[:nl]

        sig=sig.lincombo(self.cleaning)
        calib(sig) #apply calibration to cleaned PS
            
        return sig

    def get_cl_model(self,p,
                     model=None,
                     get_cmb=True,
                     get_egfs=True,
                     get_gal=True,
                     get_subpix=True,
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
                if get_gal:
                    amp = p.get(('mspec','gal','amp_%s_%s'%tuple(sorted([fr1,fr2]))),0)
                    norm_ell = float(p.get(('mspec','gal','norm_ell'),500))
                    tilt = p.get(('mspec','gal','tilt'),-1)
                    plateau_fac = p.get(('mspec','gal','plateau_fac'),1)
                    plateau_ell = p.get(('mspec','gal','plateau_ell'),300)
                    cl += amp*hstack([((plateau_ell/norm_ell)**tilt)*plateau_fac*ones(plateau_ell),(arange(plateau_ell,self.lmax)/norm_ell)**tilt])
                if get_subpix:
                    subpix = p.get(('mspec','subpix','%s_%s'%tuple(sorted([fr1,fr2]))))
                    if subpix is not None:
                        cl += subpix['amp'] * subpix['dl_shape'][:self.lmax] / subpix['dl_shape'][subpix['norm_ell']]
                    
            return model_sig

        if self.per_freq_egfs:
            model_sig = add_in_components(model_sig, in_spectra, get_cmb=get_cmb, default_egfs=False)
            model_sig = self.process_signal(p,model_sig,do_calib=False, do_beam=False)
            model_sig = add_in_components(model_sig, out_spectra, get_cmb=False, default_egfs=False)
            return model_sig.binned(self.mp['binning'])
        else:
            return self.process_signal(p,add_in_components(model_sig, in_spectra, get_cmb=get_cmb).binned(self.mp['binning']),do_calib=False, do_beam=False)


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

                   
    def get_beam_correction(self,p):
        bc = {}
        for k,v in p['mspec','beam'].items():
            ps = tuple(k.split('_'))
            bc[ps] = exp(dot(self.beampca[ps],[v.get('pca%.2i'%i,0) for i in range(self.beampca[ps].shape[1])]))
        return bc
        
    def beam_lnl(self,p):
        if self.beamcov:
            pca_vec = array([p[('mspec','beam')+tuple(k.split('.'))] for k in self.beamcov[0]])
            return dot(pca_vec,cho_solve(self.beamcov[1],pca_vec))/2
        else:
            return sum(v.get('pca%.2i'%i,0)**2/2. 
                       for k,v in p['mspec'].get('beam',{}).items()
                       for i in range(self.beampca[tuple(k.split('_'))].shape[1]))
    
    
    def lnl(self,p,model):
        
        cl_model = self.get_cl_model(p, model)
        cl_model_matrix = cl_model.get_as_matrix(lrange=self.lrange).spec[:,1]
        signal_matrix_spec = self.process_signal(p, self.signal, keep_cov=False).get_as_matrix(lrange=self.lrange)[0][:,1]
        
        #Compute the likelihood  
        dcl = cl_model_matrix - signal_matrix_spec
        return dot(dcl,cho_solve(self.signal_matrix_cov,dcl))/2 + self.beam_lnl(p)