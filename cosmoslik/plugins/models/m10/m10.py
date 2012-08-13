from cosmoslik.plugins import Model
from numpy import arange, ones, zeros

class m10(Model):
    
    def init(self, p):  
        self.ksz_template = ones(10000)
        self.tsz_template = ones(10000)
        self.clustered_template = ones(10000)
        
    def get(self, p, required):
        
        pm10 = p.get('m10',{})
        
        def get_m10_egfs(spectra, fluxcut, freqs, lmax):
            if spectra != 'cl_TT': return zeros(lmax)
            
            fr1, fr2 = freqs
            
            return sum([pm10['dgpo','amp'] * (arange(lmax)/3000.)**2 * (fr1*fr2/pm10['dgpo','norm_fr']**2)**pm10['dgpo','alpha'],
                        pm10['radio','amp'] * (arange(lmax)/3000.)**2 * (fr1*fr2/pm10['radio','norm_fr']**2)**pm10['radio','alpha'],
                        pm10['dgcl','amp'] * self.clustered_template[:lmax] * (arange(lmax)/3000.)**pm10['dgcl','tilt'] * (fr1*fr2/pm10['dgcl','norm_fr']**2)**pm10['dgcl','alpha'],
                        pm10['tsz','amp'] * self.tsz_template[:lmax],
                        pm10['ksz','amp'] * self.ksz_template[:lmax]])
    
        return {'egfs':get_m10_egfs}