import pypico
from cosmoslik_plugins.models.camb import camb
from cosmoslik import SlikPlugin

class pico(SlikPlugin):
    
    name_mapping = camb.name_mapping
    
    def __init__(self,
                 datafile):
        
        super(pico,self).__init__()
        self.pico = pypico.load_pico(datafile)
            
    def __call__(self,
                 outputs=[],
                 force=False,
                 onfail=None,
                 **kwargs):
        
        for k in list(kwargs.keys()):
            if k in self.name_mapping: 
                kwargs[self.name_mapping[k]]=kwargs[k]

        try:
            return self.pico.get(outputs=outputs,force=force,**kwargs)
        except pypico.CantUsePICO: 
            if onfail is None: raise
            else: return onfail(outputs=outputs,**kwargs)
            
