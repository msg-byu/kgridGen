! Module kpointgeneration defined in file ../src/generateKpoints.f90

subroutine f90wrap_mapkptsintobz(r, kplist, eps_, n0, n1)
    use kpointgeneration, only: mapkptsintobz
    implicit none
    
    real(8), intent(in), dimension(3,3) :: r
    real(8), intent(inout), dimension(n0,n1) :: kplist
    real(8), intent(in), optional :: eps_
    integer :: n0
    !f2py intent(hide), depend(kplist) :: n0 = shape(kplist,0)
    integer :: n1
    !f2py intent(hide), depend(kplist) :: n1 = shape(kplist,1)
    call mapkptsintobz(R=r, KpList=kplist, eps_=eps_)
end subroutine f90wrap_mapkptsintobz

subroutine f90wrap_findqpointsinzone(avecs, bvecs, n, qpoints, eps_, n0, n1)
    use kpointgeneration, only: findqpointsinzone
    implicit none
    
    real(8), intent(in), dimension(3,3) :: avecs
    real(8), intent(in), dimension(3,3) :: bvecs
    integer, intent(in) :: n
    real(8), intent(inout), dimension(n0,n1) :: qpoints
    real(8), intent(in), optional :: eps_
    integer :: n0
    !f2py intent(hide), depend(qpoints) :: n0 = shape(qpoints,0)
    integer :: n1
    !f2py intent(hide), depend(qpoints) :: n1 = shape(qpoints,1)
    call findqpointsinzone(Avecs=avecs, Bvecs=bvecs, n=n, qPoints=qpoints, &
        eps_=eps_)
end subroutine f90wrap_findqpointsinzone

! End of module kpointgeneration defined in file ../src/generateKpoints.f90

