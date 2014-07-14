#cython: cdivision=True
#cython: embedsignature=True
#cython: infertypes=True


from __future__ import division

from cosmoslik import SlikPlugin
from numpy import nan, inf
from libc.math cimport sqrt, exp, M_PI as pi, sin, sinh
from scipy.integrate import quad



cdef double taufactor = 2.30952e-5  #combination of constants (Thomson cross section, G, baryon mass) important for optical depth calc (units=Mpc^{-1})
cdef double rhoxOverOmegaxh2 = 8.09906e-11 #eV^4
cdef double G = 1.63995e2  #G in eV^(-4)Mpc^(-2)
cdef double EightPiGOver3= 8*pi/3*G  #in eV^(-4)Mpc^(-2)
#For neutrinos and photons
cdef double KelvinToeV=8.6173324e-5
cdef double Tgamma0=2.7255*KelvinToeV  
cdef double rhogamma0=pi**2/15*Tgamma0**4
cdef double OneHundredKilometersPerSecPerMpcOverSpeedofLightTimesMpc = 3.33565e-4


cdef class _cosmo_derived:

    cdef public double ombh2, omch2, omk, mnu, Nnu_massive, Nnu_massless, H0, Yp
    cdef public double ommh2, omvh2, omkh2
    
    def set_params(self, H0=None, ombh2=None, omch2=None, omk=None, mnu=None,
                   Nnu_massive=None, Nnu_massless=None, 
                   Yp=None, **kwargs):
        """
        Args:
            H0 : hubble constant today in km/s/Mpc
            ombh2, omch2, omk : densities today
            mnu : neutrino mass sum in eV
            Nnu_massive : number of massive species (mnu divided equally among them)
            Nnu_massless : number of massless species
        """
        for k,v in locals().items():
            if hasattr(self,k) and v is not None: 
                setattr(self,k,v)

        if H0 is None:
            self.omkh2 = self.ommh2 = self.omvh2 = nan
        else:
            h2 = (H0/100)**2
            self.omkh2 = omk*h2
            self.ommh2 = self.ombh2 + self.omch2
            self.omvh2 = h2 - self.omkh2 - self.ommh2 
        
    
    cpdef double Hubble(self, double z):
        """
        Returns : hubble rate at redshift z in 1/Mpc
        """
        cdef double ae, a, Te, Nphotoneff, rhonu
        
        ae=1e-5
        Te=Tgamma0/ae*(4/11)**(1/3) #scale photon temp, then convert to neutrino temp
        Nphotoneff=1+self.Nnu_massless*7/8*(4/11)**(4/3) #effective number of photon species for massless radiation
        a = 1/(1+z)
        
        if self.mnu==0: rhonu = rhogamma0*7/8*(4/11)**(4/3)*a**-4.
        else: rhonu = self.rhoredshift(a,ae,Te,self.mnu/self.Nnu_massive)
            
        return sqrt(EightPiGOver3*(rhoxOverOmegaxh2*(self.ommh2*a**-3 + self.omkh2*a**-2 + self.omvh2) 
                                   + Nphotoneff*rhogamma0*a**-4 + self.Nnu_massive*rhonu ))


    def rhoredshift(self, double a, double ae, double Te, double m):
        """
        Args:
            ae : early scale factor when still relativistic,
            Te : temperature at that time in eV
            m : mass in eV
            a : scale factor at epoch for which rho is desired
        
        Returns: rho (in units of ev^4)
        """
        
        def integrand(double x):
            return (x*x)*sqrt((x*x)*(ae/a)*(ae/a)+(m*m))/(exp(sqrt((x*x)+(m*m))/Te)+1)

        return 1/(2*pi**2)*(ae/a)**3*quad(integrand,0,inf)[0]


    def r_s(self, double z):
        """
        Returns: comoving sound horizon scale at redshift z in Mpc.
        """
        cdef double Rovera=3*self.ombh2*rhoxOverOmegaxh2/(4*rhogamma0)
        return quad(lambda double zp: 1/self.Hubble(zp) / sqrt(3*(1+Rovera/(1+zp))),z,inf)[0]

    
    cpdef double D_A(self, double z):
        """
        Returns : comoving angular-diameter distance to redshift z in Mpc.
        """
        cdef double K, dist
        K=-self.omkh2*(OneHundredKilometersPerSecPerMpcOverSpeedofLightTimesMpc)**2
        dist = self.D_prop(z)
        if K<0: return 1./sqrt(-K)*sin(dist*sqrt(-K))
        elif K > 0: return 1./sqrt(K)*sinh(dist*sqrt(K))
        elif K==0: return dist


    def D_prop(self, double z):
        """
        Returns : proper comoving distance to redshift z in Mpc.
        """
        return quad(lambda double zp: 1/self.Hubble(zp),0,z)[0]

    
    cpdef double Dv(self, double z):
        return (self.D_A(z)**2*z/self.Hubble(z))**(1/3)


#     def taub(self):
#         """ The integrand for the taub integration. """

#         cdef double Rovera=3*ombh2*rhoxOverOmegaxh2/(4*rhogamma0)  #where R = 3\rho_b/(4\rho_\gamma)

#         @vectorize
#         def taub(z,zstart=0):
#             return taufactor*(1.-Yp) * quad(lambda zp: Xe(zp)*(1.+zp)**2/(1.+Rovera/(1+zp))/Hubble(z=zp,allkw=allkw,**allkw),zstart,z,full_output=1)[0]

#         return taub

    
    
class cosmo_derived(SlikPlugin):
    
    def __init__(self):
        super(SlikPlugin,self).__init__()
        self._cosmo_derived = _cosmo_derived()
        
        for k in dir(self._cosmo_derived):
            if not k.startswith('_'):
                v = getattr(self._cosmo_derived,k)
                if callable(v): self[k] = v
