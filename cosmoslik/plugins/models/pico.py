import pypico
import cosmoslik.plugins.models.camb as camb
from cosmoslik.plugins import Model

class pico(Model):
    
    num_pico = 0
    num_camb = 0
    
    def init(self,p):
        try: datafile = p['pico','datafile']
        except KeyError: raise Exception("Please specify [pico]{datafile = ... }")
        else: self.pico = pypico.load_pico(*([datafile] if isinstance(datafile,str) else datafile))
            
        self.camb = camb.camb()
        self.camb.init(p)
        
    def get(self,p,required):
        if (self.num_pico+self.num_camb)%10==0 and p.get('pico_verbose',False): print 'PICO=%i CAMB=%i'%(self.num_pico,self.num_camb)
        try: 
            r = self.pico.get(outputs=[r for r in required if r in self.pico.outputs()],**p)
            self.camb.get(p,['z_drag'])
            self.num_pico+=1
            return r
        except pypico.CantUsePICO as e: 
            self.num_camb+=1
            print e
            print 'Calling CAMB...'
            return self.camb.get(p,required)
