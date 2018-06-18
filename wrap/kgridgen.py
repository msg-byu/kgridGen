import _kgridgen
import f90wrap.runtime
import logging

class Kpointgeneration(f90wrap.runtime.FortranModule):
    """
    Module kpointgeneration
    
    
    Defined at ../src/kpointgeneration.f90 lines 1-570
    
    """
    @staticmethod
    def mapkptsintobz(r, kplist, eps_=None):
        """
        mapkptsintobz(r, kplist[, eps_])
        
        
        Defined at ../src/kpointgeneration.f90 lines 20-105
        
        Parameters
        ----------
        r : float array
        kplist : float array
        eps_ : float
        
        """
        _kgridgen.f90wrap_mapkptsintobz(r=r, kplist=kplist, eps_=eps_)
    
    @staticmethod
    def findqpointsinzone(avecs, bvecs, n, qpoints, eps_=None):
        """
        findqpointsinzone(avecs, bvecs, n, qpoints[, eps_])
        
        
        Defined at ../src/kpointgeneration.f90 lines 541-570
        
        Parameters
        ----------
        avecs : float array
        bvecs : float array
        n : int
        qpoints : float array
        eps_ : float
        
        """
        _kgridgen.f90wrap_findqpointsinzone(avecs=avecs, bvecs=bvecs, n=n, \
            qpoints=qpoints, eps_=eps_)
    
    _dt_array_initialisers = []
    

kpointgeneration = Kpointgeneration()

