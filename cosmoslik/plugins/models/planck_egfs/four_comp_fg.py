from cosmoslik.plugins import Model
from .. egfs import egfs
from numpy import arange, loadtxt, hstack, pi, exp, zeros, ndarray
import os

class four_comp_fg(egfs):
    """
    
    """
    
    def init(self, p): 
        
        self.dir = os.path.dirname(os.path.abspath(__file__))
        padding = zeros(10000)
        def todl(cl,norm=3000): return (lambda dl: dl/dl[norm])((lambda l: cl*l*(l+1)/2/pi)(arange(cl.shape[0])))
        self.clustered_template =  todl(hstack([[0,0],loadtxt(os.path.join(self.dir,"clustered_150.dat"))[:,1],padding]))
        self.tsz_template = todl(hstack([[0,0],loadtxt(os.path.join(self.dir,"tsz.dat"))[:,1],padding]))
        
    def get_colors(self, p):
        return {'ps':'g','cl_lin':'g','cl_nonlin':'orange','tsz':'magenta'}

    def get_egfs(self, p, spectra, fluxcut, freqs, lmax, **kwargs):
        if spectra != 'cl_TT': return zeros(lmax)
        p_egfs = p.get('egfs',{})
        comps = {'ps':          p_egfs['Aps'] * (arange(lmax)/3000.)**2,
                 'cl_lin':      p_egfs['Acl_lin'] * self.clustered_template[:lmax],
                 'cl_nonlin':   p_egfs['Acl_nonlin'] * (arange(lmax)/3000.)**p_egfs.get('tilt',0.8),
                 'tsz':         p_egfs['Atsz'] * self.tsz_template[:lmax]}
        
        return comps
