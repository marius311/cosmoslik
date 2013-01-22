from cosmoslik.plugins import Model
import os
from tempfile import mkdtemp
from numpy import interp, log, exp, genfromtxt, vectorize

class recfast(Model):
    
    def init(self, p):
        self.workdir = os.path.abspath(p['recfast','workdir'] if ('recfast','workdir') in p else mkdtemp())
        self.resultfile = os.path.join(self.workdir,'result')
        self.outputfile = os.path.join(self.workdir,'output')
        self.recfast = os.path.join(os.path.abspath(os.path.dirname(__file__)),'recfast')
    
    def get(self, p, required):
           
        if 'Xe' in required:
            try: os.remove(self.resultfile)
            except: pass
            execstr = '(echo "%s \n %s %s %s \n %s %s %s \n %s \n %s" | %s) > %s' % \
                (self.resultfile,
                 p['omb'],p['omc'],p['omv'],
                 p['H0'],p['Tcmb'],p['Yp'],
                 1,
                 6,
                 self.recfast,
                 self.outputfile)
            os.system(execstr)

            try: 
                dat = genfromtxt(self.resultfile)[-2:1:-1]
                return {'Xe': lambda z: exp(interp(log(z),*log(dat.T)))}
            except:
                try: output = '\n'.join([l for l in open(self.outputfile)])
                except: output = ''
                raise Exception("RECFAST error."+output)
                
        else:
            return {}


recfast_instance = None

def xe_z(omb,omc,omv,H0,Tcmb,Yp,**kw):
    global recfast_instance
    r = recfast_instance = recfast_instance or recfast()
    r.init({})
    return r.get({'omb':omb,'omc':omc,'omv':omv,'H0':H0,'Tcmb':Tcmb,'Yp':Yp},['Xe'])['Xe']
    
    