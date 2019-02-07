MODULE wrap_kpoints
  use num_types
  use kpointgeneration, only: generateFullKpointList, generateIrredKpointList
  implicit none
CONTAINS

  !!<summary> Takes two lattices, a generating lattice (K) and the reciprocal lattice (R),
  !! and generates all the points of K that are inside one unit cell of R. </summary>
  !!<parameter name="K" regular="true"> 3x3 matrix. Columns are the generating vectors of
  !! the k-grid </parameter>
  !!<parameter name="R" regular="true"> 3x3 matrix. Columns are the reciprocal lattice
  !! vectors </parameter>
  !!<parameter name="kLVshift" regular="true"> Fractional shift of the k-grid. Given as
  !! fractions of the generating k-grid vectors. </parameter> 
  !!<parameter name="KpList" regular="true"> List of unreduced k-points in Cartesian
  !! coordinates. </parameter>
  SUBROUTINE getFullKpointList(K, R, kLVshift, KpList, reps_, aeps_)
    real(dp), intent(in) :: K(3,3)
    real(dp), intent(in) :: R(3,3)
    real(dp), intent(in) :: kLVshift(3)
    real(dp), intent(out)    :: KpList(:,:)
    real(dp), intent(in), optional:: reps_, aeps_

    real(dp), pointer :: temp_KpList(:,:)
    real(dp) :: aeps, reps

    if(.not. present(aeps_)) then
       aeps = 1e-10_dp
    else
       aeps =  aeps_
    endif

    if(.not. present(reps_)) then
       reps = 1e-10_dp
    else
       reps =  reps_
    endif
    
    KpList = 0.0_dp

    call generateFullKpointList(K, R, kLVshift, temp_KpList, reps_=reps, aeps_=aeps)

    KpList = temp_KpList
    
  end SUBROUTINE getFullKpointList


  !!<summary>Generates the reduced set of k-points and weights.</summary>
  !!<parameter name="A" regular="true">The real space lattice vectors.</parameter>
  !!<parameter name="AtBas" regular="true">The atomic basis vectors.</parameter>
  !!<parameter name="at" regular="true">The atomic site occupations.</parameter>
  !!<parameter name="K" regular="true">The k-point grid generating vectors.</parameter>
  !!<parameter name="R" regular="true">The reciprocal cell vectors.</parameter>
  !!<parameter name="kLVshift" regular="true">The k-grid shift.</parameter>
  !!<parameter name="IrrKpList" regular="true">The irreducible k-point list.</parameter>
  !!<parameter name="weights" regular="true">The weights of each of the k-points</parameter>
  !!<parameter name="reps" regular="true">Relative floating point tollerance.</parameter>
  !!<parameter name="reps" regular="true">Absolute floating point tollerance.</parameter>
  SUBROUTINE getIrredKpoints(A, AtBas, at, K, R ,kLVshift, IrrKpList, weights, reps_, aeps_)
    real(dp), intent(in) :: K(3,3), A(3,3)
    real(dp), intent(in) :: R(3,3), kLVshift(3)
    real(dp), intent(in) :: AtBas(:,:)
    integer, intent(inout) :: at(:)
    real(dp), intent(out)    :: IrrKpList(:,:)
    integer, intent(out)     :: weights(:)
    real(dp), intent(in), optional:: reps_, aeps_

    real(dp), pointer :: temp_IrrKpList(:,:)
    real(dp), allocatable :: temp_AtBas(:,:)
    integer, pointer :: temp_weights(:)
    real(dp) :: aeps, reps

    if(.not. present(aeps_)) then
       aeps = 1e-10_dp
    else
       aeps =  aeps_
    endif

    if(.not. present(reps_)) then
       reps = 1e-10_dp
    else
       reps =  reps_
    endif

    allocate(temp_AtBas(size(AtBas,1), size(AtBas,2)))
    temp_AtBas = AtBas
    
    call generateIrredKpointList(A, temp_AtBas, at, K, R, kLVshift, temp_IrrKpList, temp_weights, reps_=reps, aeps_=aeps)
    
    IrrKpList(1:size(temp_IrrKpList,1),:) = temp_IrrKpList
    weights(1:size(temp_weights,1)) = temp_weights

  end SUBROUTINE getIrredKpoints
  
end MODULE wrap_kpoints
