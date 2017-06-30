!!<summary>Auto-generated unit test for kpointgeneration.symmetryReduceKpointList
!!using FORTPY. Generated on 2017-06-29 16:39:53.598170.
!!Unit tests for          simple cubic          Bravais lattice k-points.</summary>
PROGRAM UNITTEST_symmetryReduceKpointList
  use kpointgeneration
  use num_types, only: dp
  use fortpy
  implicit none

  real(dp), pointer :: ReducedList(:,:)
  real(dp) :: shift(3)
  real(dp) :: K(3, 3)
  real(dp), allocatable :: SymOps(:,:,:)
  real(dp) :: eps_
  real(dp) :: R(3, 3)
  integer, pointer :: weights(:)
  real(dp), allocatable :: UnreducedKpList(:,:)

  real(fdp) :: fpy_start, fpy_end, fpy_elapsed = 0

  call fpy_read_f('K.in', '#', K)
  call fpy_read_f('R.in', '#', R)
  call fpy_read_f('shift.in', '#', shift)
  call fpy_read('unreduced_klist.in', '#', UnreducedKpList)
  call fpy_read('symops.in', '#', SymOps)

  call cpu_time(fpy_start)
  call symmetryReduceKpointList(K, R, shift, UnreducedKpList, SymOps, ReducedList, &
                                 weights, eps_=eps_)
  call cpu_time(fpy_end)
  fpy_elapsed = fpy_elapsed + fpy_end - fpy_start
  call pysave(weights, 'simple_cubic_wts.out')

  call pysave(fpy_elapsed, 'fpytiming.out')
END PROGRAM UNITTEST_symmetryReduceKpointList