from scipy.integrate import quad 
from numpy import *
from numpy.linalg import inv
from cosmoslik.plugins import Likelihood


class bao(Likelihood):
    """
    
    ===
    BAO
    ===
    
    Beutler et al. (2011) give rs/Dv = 0.336 \pm 0.015  at z=0.106
    Padmanabhan et al. (2012) give Dv/rs = 8.88 \pm 0.17 at z=0.35
    BOSS (Anderson et al. 2012) report Dv/rs = 13.67 \pm 0.22 at z = 0.57
    Wigglez report rs/Dv = 0.0726 \pm 0.0034 at z = 0.6
    """

    def init(self,p):
        
        #The BAO likelihoods report their constraint in terms of Dv/rs(z_drag) where
        #z_drag was calculated via the Eisenstein & Hu 1998 fitting formula, whereas
        #we actually calculate z_drag directly. This gives the correction factor to 
        #the resulting rs so that we're consistent. 
        rs_rescale = 154.6588/150.8192 
        
        self.zvec=array([0.106,0.35,0.57,0.6])
        self.Dv_over_rs_data=array([2.98,8.88,13.67,13.77])*rs_rescale
        self.Dv_over_rs_cov=(diag([0.13,0.17,0.22,0.65])*rs_rescale)**2      


    def lnl(self,p,model):
        
        rs = p['rs_drag']

        Dv_over_rs_model = [Dv(z=z, allkw=p, **p)/rs for z in self.zvec]
        
        if p.get('diagnostic'):
            from matplotlib.pyplot import ion, errorbar, plot, draw, cla, yscale, ylim, xlim
            ion()
            cla()
            errorbar(self.zvec,self.Dv_over_rs_data,yerr=sqrt(diag(self.Dv_over_rs_cov)),fmt='.')
            plot(self.zvec,Dv_over_rs_model)
            xlim(0,.7)
            draw()
        
        dx = self.Dv_over_rs_data - Dv_over_rs_model
        return dot(dx,dot(inv(self.Dv_over_rs_cov),dx))/2
        

def Dv(z, ombh2, omch2, H0, allkw, **kwargs):
    return (D_A(z=z, allkw=allkw, **allkw)**2*z/Hubble(z=z, allkw=allkw, **allkw))**(1/3.)


"""
TODO: These belong in some kind of utility module
"""

def kw_vectorize(func):
    """Like numpy.vectorize but allows keyword arguments."""
    import inspect
    vfunc = vectorize(func)
    def myfunc(**kwargs):
        return vfunc(*[kwargs[k] for k in inspect.getargspec(func).args],**kwargs)
    return myfunc




"""
TODO: All the things below here belong in a Model 
"""

taufactor = 2.30952e-5  #combination of constants (Thomson cross section, G, baryon mass) important for optical depth calc (units=Mpc^{-1})
rhoxOverOmegaxh2 = 8.09906e-11 #eV^4
G = 1.63995e2  #G in eV^(-4)Mpc^(-2)
EightPiGOver3= 8*pi/3.*G  #in eV^(-4)Mpc^(-2)
#For neutrinos and photons
KelvinToeV=8.6173324e-5
Tgamma0=2.7255*KelvinToeV  
rhogamma0=pi**2/15.*Tgamma0**4
OneHundredKilometersPerSecPerMpcOverSpeedofLightTimesMpc = 3.33565e-4


def r_s(z, ombh2, allkw, **kwargs):
    Rovera=3.*ombh2*rhoxOverOmegaxh2/(4.*rhogamma0)
    return quad(lambda zp: 1/Hubble(z=zp, allkw=allkw, **allkw) / sqrt(3*(1+Rovera/(1+zp))),z,inf)[0]
    
def D_A(z, ommh2, omkh2, omvh2, mnu, Nnu_massive, Nnu_massless, allkw, **kwargs):
    """
    dist is comoving proper distance (calculated by distance prop)
    returns comoving angular-diameter distance in Mpc
    """
    K=-omkh2*(OneHundredKilometersPerSecPerMpcOverSpeedofLightTimesMpc)**2
    dist = D_prop(z=z, allkw=allkw,**allkw)
    if K<0: return 1./sqrt(-K)*sin(dist*sqrt(-K))
    elif K > 0: return 1./sqrt(K)*sinh(dist*sqrt(K))
    elif K==0: return dist


def D_prop(z, ommh2, omkh2, omvh2, mnu, Nnu_massive, Nnu_massless, allkw, **kw):
    """returns proper comoving distance to redshift z in Mpc"""
    return quad(lambda zp: 1/Hubble(z=zp,allkw=allkw,**allkw),0,z)[0]


def Hubble(z, ommh2, omkh2, omvh2, mnu, Nnu_massive, Nnu_massless,**kw):
    """
    Returns H in Mpc^{-1}
    m is mass of single neutrino species in eV
    Nmass is number of massive species (assumed to each have same mass)
    Neff is number of massless species
    """
    ae=1.e-5
    Te=Tgamma0/ae*(4./11.)**(1./3.) #scale photon temp, then convert to neutrino temp
    Nphotoneff=1.+Nnu_massless*7./8.*(4./11.)**(4./3.) #effective number of photon species for massless radiation
    
    a = 1./(1.+z)
    
    if mnu==0: rhonu = rhogamma0*7./8.*(4./11.)**(4./3.)*a**-4.
    else: rhoredshift(a,ae,Te,mnu/(1.*Nnu_massive))

    return sqrt(EightPiGOver3*( rhoxOverOmegaxh2*( ommh2*a**(-3) + omkh2*a**(-2) + omvh2 ) + Nphotoneff*rhogamma0*a**(-4)+ Nnu_massive*rhonu ))
    
    
def rhoredshift(a, ae, Te, m):
    """
    ae is early scale factor when still relativistic,
    Te is temperature at that time in eV
    m  is mass in eV
    a  is scale factor at epoch for which rho is desired
       output:  rho (in units of ev^4)
    """
#    if ae/a*Te < 0.01*m: rho=1./(8.*pi**(3./2.))*exp(-m/Te)*m*(2*m*Te)**(1.5)*(ae/a)**3
    integral = quad(lambda x: x**2*sqrt(x**2*(ae/a)**2+m**2)/(exp(sqrt(x**2+m**2)/Te)+1.0),0.0, inf)[0]
    return 1./(2.*pi**2)*(ae/a)**3*integral


def taub(ombh2, omch2, Yp, ommh2, omkh2, omvh2, mnu, Nnu_massive, Nnu_massless, Xe, allkw, **kw):
    """ The integrand for the taub integration. """
    
    Rovera=3.*ombh2*rhoxOverOmegaxh2/(4*rhogamma0)  #where R = 3\rho_b/(4\rho_\gamma)
    
    @vectorize
    def taub(z,zstart=0):
        return taufactor*(1.-Yp) * quad(lambda zp: Xe(zp)*(1.+zp)**2/(1.+Rovera/(1+zp))/Hubble(z=zp,allkw=allkw,**allkw),zstart,z,full_output=1)[0]
    
    return taub

                   
