#cython: cdivision=True
#cython: embedsignature=True
#cython: infertypes=True
#cython: profile=False


from cosmoslik import SlikPlugin
from libc.math cimport sqrt, exp, M_PI as pi, sin, sinh, INFINITY as inf, NAN as nan
from scipy.integrate import quad, romberg, quadrature
from cosmoslik_plugins.utils.cyquad.cyquad cimport cyquad, cyquadfunc



cdef double taufactor = 2.30952e-5  #combination of constants (Thomson cross section, G, baryon mass) important for optical depth calc (units=Mpc^{-1})
cdef double rhoxOverOmegaxh2 = 8.09906e-11 #eV^4
cdef double G = 1.63995e2  #G in eV^(-4)Mpc^(-2)
cdef double EightPiGOver3= 8*pi/3*G  #in eV^(-4)Mpc^(-2)
cdef double KelvinToeV=8.6173324e-5
cdef double OneHundredKilometersPerSecPerMpcOverSpeedofLightTimesMpc = 3.33565e-4


cdef class _cosmo_derived:

    cdef public double ombh2, omch2, omk, mnu, Nnu_massive, Nnu_massless, H0, Yp, Tcmb
    cdef public double ommh2, omvh2, omkh2
    cdef public double Tgamma0, rhogamma0
    cdef public double epsrel

    def __init__(self,epsrel=1e-4):
        self.epsrel = epsrel

    def set_params(self, H0=None, ombh2=None, omch2=None, omk=None, mnu=None,
                   Nnu_massive=None, Nnu_massless=None, Tcmb=2.7255,
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

        self.Tgamma0=self.Tcmb*KelvinToeV  
        self.rhogamma0=pi**2/15*self.Tgamma0**4

        if H0 is None:
            self.omkh2 = self.ommh2 = self.omvh2 = nan
        else:
            h2 = (H0/100)**2
            self.omkh2 = omk*h2
            self.ommh2 = self.ombh2 + self.omch2
            self.omvh2 = h2 - self.omkh2 - self.ommh2 - self.rho_gamma_nu(0)/rhoxOverOmegaxh2


    cdef double rho_gamma_nu(self, double z):
        """
        Returns: energy density in neutrinos + photons in eV^4
        """
        cdef double a, Nphotoneff, rhonu
        Nphotoneff=1+self.Nnu_massless*7./8*(4./11)**(4./3) #effective number of photon species for massless radiation
        a = 1/(1+z)
        if self.mnu==0: rhonu = self.rhogamma0*7./8*(4./11)**(4./3)*a**-4
        else: rhonu = self.rhoredshift(a,self.mnu/self.Nnu_massive)
        return Nphotoneff*self.rhogamma0*a**-4 + self.Nnu_massive*rhonu
    
    cpdef double Hubble(self, double z):
        """
        Returns : hubble rate at redshift z in 1/Mpc
        """
        return sqrt(EightPiGOver3*(rhoxOverOmegaxh2*(self.ommh2*(1+z)**3 + self.omkh2*(1+z)**2 + self.omvh2) + self.rho_gamma_nu(z)))


    cdef double rhoredshift(self, double a, double m):
        """
        Args:
            ae : early scale factor when still relativistic,
            Te : temperature at that time in eV
            m : mass in eV
            a : scale factor at epoch for which rho is desired
        
        Returns: energy density in a thermal species with mass m in units of ev^4
        """
        cdef double ae=1e-5
        cdef double Te=self.Tgamma0/ae*(4./11)**(1./3) #scale photon temp, then convert to neutrino temp

        #this intergral is on an inner loop, so we do it using cyquad, which 
        # requires no Python function call overhead and is much faster. 
        return 2/(2*pi**2)*(ae/a)**3*cyquad(<cyquadfunc>rhoredshift_integrand,0,inf,self.epsrel,[a,m,ae,Te],4)

    def theta_s(self, double z):
        """
        Returns: the angular size of the sound horizon at redshift z
        """
        return self.r_s(z) / self.D_A(z)

    def r_s(self, double z):
        """
        Returns: comoving sound horizon scale at redshift z in Mpc.
        """
        cdef double Rovera=3*self.ombh2*rhoxOverOmegaxh2/(4*self.rhogamma0)
        return quad(lambda double zp: 1/self.Hubble(zp) / sqrt(3*(1+Rovera/(1+zp))),z,inf,epsabs=0,epsrel=self.epsrel)[0]

    
    cpdef double D_A(self, double z):
        """
        Returns : comoving angular-diameter distance to redshift z in Mpc.
        """
        cdef double K, dist
        K=-self.omkh2*(OneHundredKilometersPerSecPerMpcOverSpeedofLightTimesMpc)**2
        dist = self.D_prop(z)
        if K<0: return 1/sqrt(-K)*sin(dist*sqrt(-K))
        elif K > 0: return 1/sqrt(K)*sinh(dist*sqrt(K))
        elif K==0: return dist


    def D_prop(self, double z):
        """
        Returns : proper comoving distance to redshift z in Mpc.
        """
        return quad(lambda double zp: 1/self.Hubble(zp),0,z,epsabs=0,epsrel=self.epsrel)[0]

    
    cpdef double Dv(self, double z):
        return (self.D_A(z)**2*z/self.Hubble(z))**(1./3)





cdef double rhoredshift_integrand(double x, double a, double m, double ae, double Te):
    """
    The integrand for the rhoredshift integral.
    """
    return (x*x)*sqrt((x*x)*(ae/a)*(ae/a)+(m*m))/(exp(sqrt((x*x)+(m*m))/Te)+1)

    
class cosmo_derived(SlikPlugin):
    
    def __init__(self,epsrel=1e-4):
        """
        Agrs:
            epsrel : relative precision of numerical integrals (default 1e-4)
        """
        super(SlikPlugin,self).__init__()
        self._cosmo_derived = _cosmo_derived(epsrel=epsrel)

        for k in dir(self._cosmo_derived):
            if not k.startswith('_'):
                v = getattr(self._cosmo_derived,k)
                if callable(v): self[k] = v
