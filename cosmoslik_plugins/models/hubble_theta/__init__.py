from cosmoslik import SlikPlugin
            
class hubble_theta(SlikPlugin):
    
    def __init__(self):
        super(hubble_theta,self).__init__()
        from .f_hubble_theta import f_hubble_theta
        self.f_hubble_theta = f_hubble_theta
    
    def hubble_to_theta(self,
                        H0,
                        ombh2,
                        omch2,
                        omnuh2,
                        w,
                        massless_neutrinos,
                        massive_neutrinos,
                        omk,
                        **kwargs):
        h2=(H0/100.)**2
        return self.f_hubble_theta.hubble2theta(H0,ombh2/h2,omch2/h2,
                                                1-omk-(ombh2+omch2+omnuh2)/h2,omnuh2/h2,
                                                w,massless_neutrinos,massive_neutrinos)
    
    def theta_to_hubble(self,
                        theta,
                        ombh2,
                        omch2,
                        omk,
                        omnuh2,
                        w,
                        massless_neutrinos,
                        massive_neutrinos,
                        **kwargs):
        
        return self.f_hubble_theta.theta2hubble(theta,ombh2,omch2,omk,omnuh2,
                                                w,massless_neutrinos,massive_neutrinos)
        
