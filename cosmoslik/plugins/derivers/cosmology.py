from numpy import exp, log
from cosmoslik.plugins import Deriver

aliases = [
     ['As','A_s','scalar_amp(1)'],
     ['ns','n_s','Ns','scalar_spectral_index(1)'],
     ['nrun','n_run','scalar_nrun(1)'],
     ['omb','omegab','omega_baryon','om_b','omega_b'],
     ['omc','omegac','omega_c'],
     ['omn','omegan','omega_n','omnu','omeganu','omega_nu'],
     ['omv','omegav','om_v','omega_v','omega_lambda'],
     ['omk','omegak','omega_k'],
     ['ombh2','omegabh2'],
     ['omch2','omegach2'],
     ['omkh2','omegakh2'],
     ['omvh2','omegavh2','omlh2','omegalh2'],
     ['omnh2','omeganh2','omnuh2','omeganuh2'],
     ['H0','H','hubble'],
     ['theta','thetastar'],
     ['Nnu_massless','massless_neutrinos','Neff'],
     ['Nnu_massive','massive_neutrinos'],
     ['Yp','helium_fraction'],
     ['tau','re_optical_depth']
]

def add_aliases(aliases,p):
    """Adds keys to dictionary corresponding to list of aliases"""
    for a in aliases:
        h = {k:p[k] for k in a if k in p}
        if any(h):
            if not all(h.values()[0]==v or abs(h.values()[0]/v-1)<1e-4 for v in h.values()): raise Exception("You have aliased keys with different values: "+str(h))
            for k in a: p[k]=h.values()[0]
            
            
class cosmology(Deriver):


    def add_derived(self,p):
        """
        Given a dictionary of cosmological parameters attempts to add as many
        aliases and derived parameters as possible
        """
        from hubble_theta import hubble2theta, theta2hubble
        add_aliases(aliases,p)
        if 'theta' in p:
            if all(k in p for k in ['ombh2','omch2', 'omk', 'omnuh2', 'w', 'Nnu_massless', 'Nnu_massive']):
                h=p["h"]=theta2hubble(p["theta"],p["ombh2"],p["omch2"],p["omk"],p["omnuh2"],p["w"],p["Nnu_massless"],p["Nnu_massive"])/100.
                for k in ['omb','omc','omm','omnu','omv','omk']:
                    if k+'h2' in p: p[k]=p[k+'h2']/h**2
                    elif k in p: p[k+'h2']=p[k]*h**2
                p["omv"]=1-p['omb']-p['omc']-p['omnu']
                p["omvh2"]=p["omv"]*h**2
                p["omm"]=p["omb"]+p["omc"]
                p["ommh2"]=p["omm"]*h**2
                p['mnu'] = p['omnuh2']/94.
                p["H0"]=100*h
        elif 'H0' in p or 'h' in p:
            if 'H0' in p: h=p['h']=p['H0']/100.
            else: p['H0']=p['h']*100; h=p['h']  
            for k in ['omb','omc','omm','omnu','omv','omk']:
                if k in p: p[k+'h2']=p[k]*h**2
                elif k+'h2' in p: p[k]=p[k+'h2']/h**2
            if 'omv' not in p: p['omv']=1-p['omb']-p['omc']-p['omnu']
            if all(k in p for k in ['omb','omc', 'omv', 'omn', 'w', 'Nnu_massless', 'Nnu_massive']):
                p['theta']=hubble2theta(p['H0'],p['omb'],p['omc'],p['omv'],p['omn'],p['w'],p['Nnu_massless'],p['Nnu_massive'])
        
        if 'As' in p: p['logA']=log(10**10*p['As'])
        elif 'logA' in p: p['As']=exp(p['logA'])*10**(-10)
        add_aliases(aliases,p)
