import _kpoints
import f90wrap.runtime
import logging
import numpy as np

def get_irr_kpoints(atoms, grid=None, HNF=None, shift=None, reps=None, aeps=None):
    """Generates the irreducible k-point set for the system defined by the
    atoms object and either the grid or the HNF given. If both the HNF
    and grid are provided then the grid will be used and the HNF ignored.

    Args:
        atoms (ase.Atoms): an ASE Atoms object.
        HNF (optional, matrix, int): an integer matrix that generates a supercell of 
          the parent cell.
        grid (optional, matrix. float): the grid generating vectors, as columns of 
          the matrix.
        shift (optional, list, float): a list of the offset from gamma.

    Returns:
        The irreducible k-points and their weights.
    """

    if shift is None:
        shift = [0,0,0]
    if reps is None:
        reps = 1E-10
    if aeps is None:
        aeps = 1E-10
    
    if HNF is None and grid is None:
        raise ValueError("Either the HNF or the grid generating vectors must "
                         "be provided.")

    cell = np.transpose(atoms.cell)
    at_base = []
    for pos in atoms.positions:
        at_base.append(np.matmul(np.linalg.inv(cell),pos))
    cell = cell/np.linalg.det(cell)
    reciprocal_cell = np.transpose(np.linalg.inv(cell))
    if grid is not None:
        int_mat = np.matmul(np.linalg.inv(grid),reciprocal_cell)
        vol = np.linalg.det(int_mat)
        if np.allclose(vol, np.round(vol)):
            vol = int(np.round(vol))
        else:
            raise ValueError("The grid and the reciprocal lattice aren't commensurate.")

    if HNF is not None and grid is None:
        vol = int(np.round(np.linalg.det(HNF)))
        supercell = np.matmul(cell,HNF)
        grid = np.linalg.inv(np.transpose(supercell))


    at_base = np.transpose(at_base)
    t_map = {k:v for v,k in enumerate(np.unique(atoms.get_atomic_numbers()))}
    at = np.array([t_map[i] for i in atoms.get_atomic_numbers()], dtype=np.int32, order='F')
    irr_kpoint_list = np.zeros((vol,3), dtype=float, order='F')
    weights = np.zeros(vol,dtype=np.int32, order='F')

    print("grid",grid)
    print("grid det", np.linalg.det(grid))
    Wrap_Kpoints.getirredkpoints(cell, at_base, at, grid, reciprocal_cell,
                                            shift, irr_kpoint_list, weights, reps_=reps, aeps_=aeps)

    keep = np.where(weights != 0)[0]

    return irr_kpoint_list[keep], weights[keep]
    
def get_all_kpoints(atoms, HNF=None, grid=None, shift=None, eps=None, aeps=None):
    """Generates the full k-points list for the system defined by the
    atoms object and either the grid or the HNF given. If both the HNF
    and grid are provided then the grid will be used and the HNF ignored.

    Args:
        atoms (ase.Atoms): an ASE Atoms object.
        HNF (optional, matrix, int): an integer matrix that generates a supercell of 
          the parent cell.
        grid (optional, matrix. float): the grid generating vectors, as columns of 
          the matrix.
        shift (optional, list, float): a list of the offset from gamma.

    Returns:
        A full list of the k-points in the reciprocal cell.
    """

    if shift is None:
        shift = [0,0,0]
    if eps is None:
        reps = 1E-10
    if aeps is None:
        aeps = 1E-10
    
    if HNF is None and grid is None:
        raise ValueError("Either the HNF or the grid generating vectors must "
                         "be provided.")
    
    reciprocal_cell = np.linalg.inv(atoms.cell)
    if grid is not None:
        int_mat = np.matmul(np.linalg.inv(grid),reciprocal_cell)
        vol = np.linalg.det(int_mat)
        if np.allclose(vol, np.round(vol)):
            vol = int(np.round(vol))
        else:
            raise ValueError("The grid and the reciprocal lattice aren't commensurate.")

    if HNF is not None and grid is None:
        vol = int(np.round(np.linalg.det(HNF)))
        supercell = np.matmul(np.transpose(atoms.cell),HNF)
        grid = np.linalg.inv(np.transpose(supercell))

    kpoint_list = np.zeros((vol,3), dtype=float, order='F')

    Wrap_Kpoints.getfullkpointlist(grid, reciprocal_cell,
                                   shift, kpoint_list, reps_=reps, aeps_=aeps)

    return kpoint_list


class Wrap_Kpoints(f90wrap.runtime.FortranModule):
    """
    Module wrap_kpoints
    
    
    Defined at ../src/wrap_kpoints.f90 lines 1-92
    
    """
    @staticmethod
    def getfullkpointlist(k, r, klvshift, kplist, reps_, aeps_):
        """
        getfullkpointlist(k, r, klvshift, kplist, reps_, aeps_)
        
        
        Defined at ../src/wrap_kpoints.f90 lines 17-46
        
        Parameters
        ----------
        k : float array
        r : float array
        klvshift : float array
        kplist : float array
        reps_ : float
        aeps_ : float
        
        """
        _kpoints.f90wrap_getfullkpointlist(k=k, r=r, klvshift=klvshift, kplist=kplist, \
            reps_=reps_, aeps_=aeps_)
    
    @staticmethod
    def getirredkpoints(a, atbas, at, k, r, klvshift, irrkplist, weights, reps_, \
        aeps_):
        """
        getirredkpoints(a, atbas, at, k, r, klvshift, irrkplist, weights, reps_, aeps_)
        
        
        Defined at ../src/wrap_kpoints.f90 lines 59-92
        
        Parameters
        ----------
        a : float array
        atbas : float array
        at : int array
        k : float array
        r : float array
        klvshift : float array
        irrkplist : float array
        weights : int array
        reps_ : float
        aeps_ : float
        
        """
        _kpoints.f90wrap_getirredkpoints(a=a, atbas=atbas, at=at, k=k, r=r, \
            klvshift=klvshift, irrkplist=irrkplist, weights=weights, reps_=reps_, aeps_=aeps_)
    
    _dt_array_initialisers = []
    

wrap_kpoints = Wrap_Kpoints()

