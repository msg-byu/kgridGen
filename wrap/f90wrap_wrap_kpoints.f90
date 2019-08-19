! Module wrap_kpoints defined in file ../src/wrap_kpoints.f90

subroutine f90wrap_getfullkpointlist(k, r, klvshift, kplist, reps_, aeps_, n0, &
    n1)
    use wrap_kpoints, only: getfullkpointlist
    implicit none
    
    real(8), intent(in), dimension(3,3) :: k
    real(8), intent(in), dimension(3,3) :: r
    real(8), intent(in), dimension(3) :: klvshift
    real(8), intent(inout), dimension(n0,n1) :: kplist
    real(8) :: reps_
    real(8) :: aeps_
    integer :: n0
    !f2py intent(hide), depend(kplist) :: n0 = shape(kplist,0)
    integer :: n1
    !f2py intent(hide), depend(kplist) :: n1 = shape(kplist,1)
    call getfullkpointlist(K=k, R=r, kLVshift=klvshift, KpList=kplist, reps_=reps_, &
        aeps_=aeps_)
end subroutine f90wrap_getfullkpointlist

subroutine f90wrap_getirredkpoints(a, atbas, at, k, r, klvshift, irrkplist, &
    weights, reps_, aeps_, n0, n1, n2, n3, n4, n5)
    use wrap_kpoints, only: getirredkpoints
    implicit none
    
    real(8), intent(in), dimension(3,3) :: a
    real(8), intent(in), dimension(n0,n1) :: atbas
    integer, intent(inout), dimension(n2) :: at
    real(8), intent(in), dimension(3,3) :: k
    real(8), intent(in), dimension(3,3) :: r
    real(8), intent(in), dimension(3) :: klvshift
    real(8), intent(inout), dimension(n3,n4) :: irrkplist
    integer, intent(inout), dimension(n5) :: weights
    real(8) :: reps_
    real(8) :: aeps_
    integer :: n0
    !f2py intent(hide), depend(atbas) :: n0 = shape(atbas,0)
    integer :: n1
    !f2py intent(hide), depend(atbas) :: n1 = shape(atbas,1)
    integer :: n2
    !f2py intent(hide), depend(at) :: n2 = shape(at,0)
    integer :: n3
    !f2py intent(hide), depend(irrkplist) :: n3 = shape(irrkplist,0)
    integer :: n4
    !f2py intent(hide), depend(irrkplist) :: n4 = shape(irrkplist,1)
    integer :: n5
    !f2py intent(hide), depend(weights) :: n5 = shape(weights,0)
    call getirredkpoints(A=a, AtBas=atbas, at=at, K=k, R=r, kLVshift=klvshift, &
        IrrKpList=irrkplist, weights=weights, reps_=reps_, aeps_=aeps_)
end subroutine f90wrap_getirredkpoints

! End of module wrap_kpoints defined in file ../src/wrap_kpoints.f90

