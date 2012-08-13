import os, sys, re
from numpy import zeros, loadtxt, hstack
from cosmoslik.plugins import Model

try: import camb4py
except ImportError: import local_camb4py as camb4py

class camb(Model):
    """
    ====
    CAMB
    ====
    
    - Written by Antony Lewis and Anthony Challinor (see `<http://camb.info>`_)
    - CosmoSlik module by Marius Millea
    
    
    Usage
    =====
    
    To use this module, add ``camb`` to the list of ``models``.
    
    Notes
    =====
    
    CosmoSlik uses a modified CAMB which under-the-hood reads the ini file 
    from `stdin` and writes results to `stdout`. This module also handles some
    of the logic in chosing settings. In particular,
    
    - It automatically choses which of ``get_scalar``, ``get_tensor``, 
      ``get_vector``, and ``get_transfer`` should be set to true based 
      on which likelihoods are used. 
    - It automatically sets ``do_lensing=True`` if ``Alens>0``.
    - If automatically increases ``lmax`` used in calculations to ensure 
      full accuracy of outputs up to the specified ``lmax``. 
    
    .. note:: 
    
        For users wishing to modify or add functionality to CAMB, we recommend starting
        by making a copy of this module and leaving ``inidriver.F90`` mostly
        intact (it contains the under-the-hood changes mentioned above.)
    
    Parameters
    ==========
    
    [camb].X
    --------
        X can be any valid CAMB parameter that you would normally give in the
        CAMB ini file. These parameters are passed straight to CAMB (and can
        be sampled over).  
        
    [camb].defaults
    ---------------
        Path to a CAMB ini file which specifies default values for parameters not set 
        by the user or by CosmoSlik. (default: defaults.ini in the module folder) 
        
    Alens
    -----
        cl = Alens * cl_lensed + (1-Alens) * cl_scalar 
    
    
    """
    
    def init(self,p):
        pcamb = p.get('camb',{})
        self.cambdir = os.path.abspath(os.path.join(os.path.dirname(__file__),'camb'))
        pcamb.setdefault('executable',os.path.join(self.cambdir,'camb'))
        if not os.path.exists(pcamb['executable']): raise Exception("Could not find camb executable at '%s'"%pcamb['executable'])
        self.cambdefs = camb4py.read_ini(pcamb.get('defaults',os.path.join(self.cambdir,'defaults.ini')))
        self.cambdefs.update(pcamb)
        self.camb = camb4py.load(pcamb['executable'], self.cambdefs)
    
    
    def get(self,p,required):
        pcamb = p.get('camb',{})
    
        cambini = pcamb['ini'] = dict(self.cambdefs)
        cambini.update(pcamb)
        cambini.update(p)
        
        Alens = p.get('Alens',1)
        cambini['get_scalar_cls'] = doscal = any(x in required for x in ['cl_TT','cl_TE','cl_EE','cl_BB','cl_pp','cl_pT'])
        cambini['get_tensor_cls'] = dotens = (p.get('r',0) != 0)
        cambini['get_transfer'] = dotrans = any(x in required for x in ['lin_mpk','nonlin_mpk','trans'])
        if 'nonlin_mpk' in required: cambini['do_nonlinear'] = min(1,cambini.get('do_nonlinear',1))
        cambini['do_lensing'] = dolens = (doscal and Alens != 0)
        docl = doscal or dolens or dotens 
        if docl:
            lmax = pcamb['lmax']
            cambini['l_max_scalar'] = lmax + 50 + (100 if dolens else 0)
            lmax_tens = cambini['l_max_tensor'] = p.get('lmax_tensor',lmax + 50)
        
        for k,v in cambini.items():
            if not isinstance(v,(float, int, str, bool)): cambini.pop(k)
        
        result = {}

        #Call CAMB
        output = self.camb(**cambini)

        try:
        
            if doscal: scal = dict(zip(['l','TT','EE','TE','pp','pT'],output['scalar'].T))
            if dolens: lens = dict(zip(['l','TT','EE','BB','TE'],output['lensed'].T))
            if dotens: tens = dict(zip(['l','TT','EE','BB','TE'],output['tensor'].T))
            if dotrans: 
                for x in ['lin_mpk','nonlin_mpk','trans']:
                    if x in required: result[x]=output[x]
                    
            #Combine cl contributions
            if docl:
                for x in ['TT','TE','EE','BB']: 
                    if 'cl_%s'%x in required:
                        result['cl_%s'%x] = zeros(lmax)
                        if doscal or dolens: 
                            result['cl_%s'%x][2:lmax] += (((1-Alens)*scal[x][:lmax-2] if x!='BB' and doscal else 0)) + (Alens*lens[x][:lmax-2] if dolens else 0)
                        if dotens:
                            result['cl_%s'%x][2:lmax_tens] += tens[x][:lmax_tens-2]
                if dolens:
                    if 'cl_pp' in required: result['cl_pp'] = hstack([[0,0],scal['pp'][:lmax-2]])
                    if 'cl_pT' in required: result['cl_pT'] = hstack([[0,0],scal['pT'][:lmax-2]])

            #TODO: figure out where to put this stuff
            p['z_drag'] = float(output['misc']['z_drag'])
            p['rs_drag'] = float(output['misc']['rs_drag'])
        
        except Exception as e:
            raise Exception("""An error occurred reading CAMB result: %s \nCAMB output:\n"""%repr(e)+output['stdout'])


        return result

