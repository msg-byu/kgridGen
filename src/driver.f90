PROGRAM kpoint_driver
  use kpointGeneration
  use num_types
  use vector_matrix_utilities
  use symmetry, only : get_lattice_pointGroup
  use rational_mathematics, only: HermiteNormalForm
  use fortpy, only : pysave
  implicit none

  real(dp)              :: K(3,3), R(3,3), Hinv(3,3), eps, shift(3), M(3,3), minkedR(3,3), minkedK(3,3), minked_shift(3)
  real(dp), pointer     :: klist(:,:)
  real(dp), pointer     :: pgOps(:,:,:), rdKlist(:,:)
  integer, pointer      :: weights(:)
  integer H(3,3), i
  
  ! Finite precision tolerance (same as default value)
  eps = 1e-6_dp
  
  ! Reciprocal lattice vectors
  ! R = transpose(reshape((/ 1.98520863_dp,    0.00000000_dp,   0.00000000_dp, &
  !                          0.0_dp,           1.44640546_dp,   0.00000000_dp, &
  !                          0.0575324872_dp,  0.0_dp,          1.42600347_dp /),(/3,3/)))

  ! HNF Matrix
  ! H = transpose(reshape((/ 2, 0, 0, &
  !                          0, 2, 0, &
  !                          0, 0, 2 /),(/3,3/)))
  
  ! Shift of grid in grid coordinates.
  ! shift = (/ 5.8_dp, 5.8_dp, 5.8_dp /)

  shift = (/ 0.5_dp, -4.5_dp, 0.5_dp /)
  R = transpose(reshape((/ 1_dp, 0_dp, 0_dp, &
                          10_dp, 1_dp, 0_dp, &
                          0_dp, 0_dp, 1_dp /),(/3,3/)))
  H = transpose(reshape((/ 4, 0, 0, &
                          0, 4, 0, &
                          0, 0, 4 /),(/3,3/)))

  call matrix_inverse(real(H,dp), Hinv, eps_=1e-12_dp)
  
  ! Columns of K are the grid generating vectors.
  K = matmul(R,Hinv)
  
  write(*,'(3("R: ",3(1x,f7.3),/))') (R(i,:),i=1,3)
  write(*,'(3("H: ",3(1x,i3),/))') (H(i,:),i=1,3)
  Write(*,'(3("Hinv: ",3(1x,f7.3),/))') (Hinv(i,:),i=1,3)
  write(*,'(3("K: ",3(1x,f7.3),/))') (K(i,:),i=1,3)
  write(*,'("shift: ",3(f6.3,1x))') shift
  write(*,'("cart shift: ",3(f6.3,1x))') matmul(K,shift)

  ! write(*,'(3("PP: ",3(1x,f7.3),/))') matmul(K,(/1,0,0/))  
  call generateFullKpointList(K, R, shift, klist)
  do i = 1,determinant(H)
     write(*,'(3(1x,g11.4))') klist(i,:)
  end do
  
  ! call get_lattice_pointGroup(R, pgOps, eps)
  ! Normal tests
  ! call pysave(K, "../tests/simple_cubic/K.in.10")
  ! call pysave(R, "../tests/simple_cubic/R.in.10")
  ! call pysave(shift, "../tests/simple_cubic/shift.in.10")
  ! call pysave(klist, "../tests/simple_cubic/klist.in.10")
  ! call pysave(pgOps, "../tests/simple_cubic/symops.in.10")
  ! call symmetryReduceKpointList(K, R, shift,  klist, pgOps, rdKlist, weights, eps)
  ! call pysave(rdKlist, "../tests/simple_cubic/simple_cubic_kpts.out.10")
  ! call pysave(weights, "../tests/simple_cubic/tetragonal_wts.out.10")

  ! VASP
  ! call pysave(K, "../tests/body-centered_tetragonal/K.in.vasp10")
  ! call pysave(R, "../tests/body-centered_tetragonal/R.in.vasp10")
  ! call pysave(shift, "../tests/body-centered_tetragonal/shift.in.vasp10")
  ! call pysave(klist, "../tests/body-centered_tetragonal/klist.in.vasp10")
  ! call pysave(pgOps, "../tests/body-centered_tetragonal/symops.in.vasp10")
  ! call symmetryReduceKpointList(K, R, shift,  klist, pgOps, rdKlist, weights, &
  !      eps_=eps)
  ! call pysave(rdKlist, "../tests/body-centered_tetragonal/klist.out.vasp10")
  ! call pysave(weights, "../tests/body-centered_tetragonal/weights.out.vasp10")
  
  ! call pysave(R, "../tests/simple_cubic/rlatvecs.in.BZM1")
  ! call pysave(klist, "../tests/simple_cubic/klist.in.BZM1")
  
  ! Map into first Brillouin zone tests.
  call mapKptsIntoMinkowskiUnitCell(R, klist, M, eps)

  minkedK = matmul(M,K)
  minkedR = matmul(M,R)
  minked_shift = matmul(M,shift)
  
  call get_lattice_pointGroup(minkedR, pgOps, eps)
  
  call symmetryReduceKpointList(minkedK, minkedR, minked_shift, klist, pgOps, &
       rdKlist, weights, eps_=.0001_dp)

  call mapKptsIntoFirstBZ(minkedR, rdKlist, eps)

  ! call symmetryReduceKpointList(K, R, shift, klist, pgOps, rdKlist, weights, &
  !      eps_=eps)

  ! call mapKptsIntoFirstBZ(R, rdKlist, eps)
  write(*,'(//"**********")')
  write(*,'(//"**********")')
  
  do i = 1,size(rdKlist,1)
     write(*,'(3(1x,f7.3),3x,"w:",i5)') rdKlist(i,:),weights(i)
  end do

  ! call pysave(klist, "../tests/simple_cubic/klist.out.BZM1")
  ! write(*,'(3("minkedR: ",3(1x,f7.3),/))') (minkedR(i,:),i=1,3)
!
!  write(*,'(//)')
!  write(*,'("Unrd kpts: ",i7)') size(klist,1)
!  write(*,'("Rdcd kpts: ",i7)') size(rdKlist,1)
!  write(*,'("Rdn ratio: ",3x,f4.1)') size(klist,1)/real(size(weights))

END PROGRAM kpoint_driver
