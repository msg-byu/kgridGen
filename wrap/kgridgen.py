import _kgridgen
import f90wrap.runtime
import logging

class Kpointgeneration(f90wrap.runtime.FortranModule):
    """
    Module kpointgeneration
    
    
    Defined at ../src/kpointgeneration.f90 lines 1-570
    
    """
    @staticmethod
    def mapkptsintobz(r, kplist, reps_=None, aeps_=None):
        """
        mapkptsintobz(r, kplist[, eps_])
        
        
        Defined at ../src/kpointgeneration.f90 lines 20-105
        
        Parameters
        ----------
        r : float array
        kplist : float array
        eps_ : float
        
        """
        _kgridgen.f90wrap_mapkptsintobz(r=r, kplist=kplist, reps_=reps_, aeps_=aeps_)
    
    @staticmethod
    def findqpointsinzone(avecs, bvecs, n, qpoints, reps_=None, aeps_=None):
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
                                            qpoints=qpoints, reps_=reps_, aeps_=aeps_)
    
    _dt_array_initialisers = []
    

kpointgeneration = Kpointgeneration()
