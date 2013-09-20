import os, sys, re
from numpy import zeros, loadtxt, hstack, arange
from ConfigParser import RawConfigParser
from StringIO import StringIO

from cosmoslik import SlikPlugin, SubprocessExtension


class camb(SlikPlugin):
    """
    """
    
    name_mapping = {'H0':'hubble',
                    'As':'scalar_amp(1)',
                    'ns':'scalar_spectral_index(1)',
                    'nrun':'scalar_nrun(1)',
                    'Yp':'helium_fraction',
                    'tau':'re_optical_depth'}
    
    def __init__(self):
        super(camb,self).__init__()
        
        self._cambdir = os.path.abspath(os.path.join(os.path.dirname(__file__),'camb'))
        self._cambdefs = read_ini(os.path.join(self._cambdir,'defaults.ini'))
        self._cambdefs['HighLExtrapTemplate'] = os.path.abspath(os.path.join(self._cambdir,'HighLExtrapTemplate_lenspotentialCls.dat'))
        self._camb = camb_f2py()
        pass
        
    
    def __call__(self,
                 ombh2,
                 omch2,
                 H0,
                 As,
                 ns,
                 tau,
                 omnuh2=0,
                 w=-1,
                 Alens=1,
                 r=0,
                 nrun=0,
                 omk=0,
                 massive_neutrinos=3.046,
                 massless_neutrinos=0.000,
                 l_max_scalar=3000,
                 l_max_tensor=3000,
                 outputs=[],
                 **kwargs):
        
        locs = locals()
        cambini = dict(self._cambdefs)
        cambini.update(locs)
        
        for k in cambini.keys():
            if k in self.name_mapping: 
                cambini[self.name_mapping[k]]=cambini.pop(k)
        
        cambini['get_scalar_cls'] = doscal = any(x in outputs for x in ['cl_TT','cl_TE','cl_EE','cl_BB','cl_pp','cl_pT'])
        cambini['get_tensor_cls'] = dotens = (r != 0)
        cambini['get_transfer'] = dotrans = any(x in outputs for x in ['lin_mpk','nonlin_mpk','trans'])
        if 'nonlin_mpk' in outputs: cambini['do_nonlinear'] = min(1,cambini.get('do_nonlinear',1))
        cambini['do_lensing'] = dolens = (doscal and Alens != 0)
        docl = doscal or dolens or dotens 
        
        for k,v in cambini.items():
            if not isinstance(v,(float, int, str, bool)): cambini.pop(k)
        
        result = {}

        #Call CAMB
        output = self._camb(**cambini)
        

        try:
        
            if doscal: scal = dict(zip(['l','TT','EE','TE','pp','pT'],output['scalar'].T))
            if dolens: lens = dict(zip(['l','TT','EE','BB','TE'],output['lensed'].T))
            if dotens: tens = dict(zip(['l','TT','EE','BB','TE'],output['tensor'].T))
            if dotrans: 
                for x in ['lin_mpk','nonlin_mpk','trans']:
                    if x in outputs: result[x]=output[x]
                    
            #Combine cl contributions
            if docl:
                for x in ['TT','TE','EE','BB']: 
                    if 'cl_%s'%x in outputs:
                        result['cl_%s'%x] = zeros(l_max_scalar)
                        if doscal or dolens: 
                            _lmax = min(l_max_scalar,(lens if dolens else scal)[x].shape[0])
                            result['cl_%s'%x][2:_lmax] += (((1-Alens)*scal[x][:_lmax-2] if x!='BB' and doscal else 0)) + (Alens*lens[x][:_lmax-2] if dolens else 0)
                        if dotens:
                            _lmax = min(l_max_tensor,tens[x].shape[0])
                            result['cl_%s'%x][2:_lmax] += tens[x][:_lmax-2]
                if dolens:
                    if 'cl_pp' in outputs: result['cl_pp'] = hstack([[0,0],scal['pp'][:l_max_scalar-2]])
                    if 'cl_pT' in outputs: result['cl_pT'] = hstack([[0,0],scal['pT'][:l_max_scalar-2]])

            #TODO: figure out where to put this stuff
            result['zdrag'] = float(output['misc']['z_drag'])
        
        except Exception as e:
            print e
            print output
            raise
            #raise Exception("""An error occurred reading CAMB result: %s \nCAMB output:\n"""%repr(e)+output.get('stdout',''))

        return result #TODO make result class


def read_ini(ini):
    """Read CAMB ini into dictionary"""
    if os.path.exists(ini): ini = open(ini).read()
    config = RawConfigParser()
    config.optionxform=str
    config.readfp(StringIO('[root]\n'+ini))
    return dict(config.items('root'))


def try_bool2str(value):
    if value is True: return 'T'
    elif value is False: return 'F'
    else: return value


class camb_f2py(object):
    
    def __init__(self):
        self.pycamb = SubprocessExtension('pycamb',globals())
        
    def __call__(self, **params):
        lines = ['%s = %s'%(k,try_bool2str(v)) for k,v in params.items()]
        max_line_len = max(len(l) for l in lines)
        sp = ''.join(sorted([l.ljust(max_line_len) for l in lines]))
        self.pycamb.run(sp,max_line_len,len(lines))

        def add_ell(cl):
            if cl is not None:
                nl = cl.shape[0]
                return hstack([arange(2,nl+2).reshape(nl,1),cl[:,0]])

        result = {}
        if self.pycamb.cl_scalar_ is not None: result['scalar'] = add_ell(self.pycamb.cl_scalar_)
        if self.pycamb.cl_lensed_ is not None: result['lensed'] = add_ell(self.pycamb.cl_lensed_)
        if self.pycamb.cl_tensor_ is not None: result['tensor'] = add_ell(self.pycamb.cl_tensor_)
        if self.pycamb.mpk_lin_ is not None: result['lin_mpk'] = self.pycamb.mpk_lin_[:,:,0,0].T
        if self.pycamb.mpk_nonlin_ is not None: result['nonlin_mpk'] = self.pycamb.mpk_nonlin_[:,:,0,0].T
        if self.pycamb.transfer_ is not None: result['trans'] = self.pycamb.transfer_[:,:,0].T
        result['misc'] = {}
        if self.pycamb.z_drag_ is not None: result['misc']['z_drag'] = float(self.pycamb.z_drag_)
        
        return result
    
