!!<summary>Auto-generated unit test for kpointgeneration.symmetryReduceKpointList
!!using FORTPY. Generated on 2017-07-07 15:52:12.023111.
!!Unit tests for simple cubic         Bravais lattice k-points.</summary>
PROGRAM UNITTEST_symmetryReduceKpointList
  use kpointgeneration
  use num_types, only: dp
  use fortpy
  implicit none

  real(dp) :: kLVshift(3)
  real(dp), pointer :: ReducedList(:,:)
  real(dp) :: K(3, 3)
  real(dp), allocatable :: SymOps(:,:,:)
  real(dp) :: eps_
  real(dp) :: R(3, 3)
  integer, pointer :: weights(:)
  real(dp), allocatable :: UnreducedKpList(:,:)

  real(fdp) :: fpy_start, fpy_end, fpy_elapsed = 0

  call fpy_read_f('symops.in', '#', SymOps)
  eps_ = 1E-10_dp
  call fpy_read_f('K.in', '#', K)
  call fpy_read_f('R.in', '#', R)
  call fpy_read_f('shift.in', '#', kLVshift)
  call fpy_read_f('unreduced_klist.in', '#', UnreducedKpList)

  call cpu_time(fpy_start)
  call symmetryReduceKpointList(K, R, kLVshift, UnreducedKpList, SymOps, ReducedList, &
                                 weights, eps_=eps_)
  call cpu_time(fpy_end)
  fpy_elapsed = fpy_elapsed + fpy_end - fpy_start
  call pysave(reducedlist, 'face-centered_cubic_kpts.out')

  call pysave(fpy_elapsed, 'fpytiming.out')
END PROGRAM UNITTEST_symmetryReduceKpointList