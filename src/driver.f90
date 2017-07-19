PROGRAM kpoint_driver
  use kpointGeneration
  use num_types
  use vector_matrix_utilities
  use symmetry, only : get_lattice_pointGroup
  use rational_mathematics, only: HermiteNormalForm
  use fortpy, only : pysave
  implicit none

  real(dp)              :: K(3,3), R(3,3), Hinv(3,3), eps, shift(3)
  real(dp), allocatable :: klist(:,:)
  real(dp), pointer     :: pgOps(:,:,:), rdKlist(:,:)
  integer, pointer      :: weights(:)
  integer H(3,3), i
  
  ! Finite precision tolerance (same as default value)
  eps = 1e-6_dp
  
  ! Reciprocal lattice vectors
  ! R = transpose(reshape((/ 1.98520863_dp,    0.00000000_dp,   0.00000000_dp, &
  !                          0.0_dp,           1.44640546_dp,   0.00000000_dp, &
  !                          0.0575324872_dp,  0.0_dp,          1.42600347_dp /),(/3,3/)))

  ! R = transpose(reshape((/     /),(/3,3/)))

  R = transpose(reshape((/  1.89138631_dp, 1.89138631_dp, 0.0_dp, &
                            -1.09199239_dp, 1.09199239_dp, 0.0_dp, &
                            0.0_dp, 0.0_dp, 1.21814372_dp /),(/3,3/)))
  
  ! R = transpose(reshape((/  0.0_dp,  0.8_dp,  0.8_dp, &
  !                           1.35_dp,  0.0_dp, 1.35_dp, &
  !                           1.7_dp,  1.7_dp,  0.0_dp  /),(/3,3/)))
  
  ! HNF Matrix
  H = transpose(reshape((/ 32, 0, 0, &
                           0, 32, 0, &
                           0, 0, 32 /),(/3,3/)))
  
  shift = (/ 0.0_dp, 0.0_dp, 0.0_dp /)
  ! shift = (/  2.0_dp/3.0_dp, 2.0_dp/3.0_dp, 2.0_dp/3.0_dp /)
  
  call matrix_inverse(real(H,dp), Hinv, eps_=1e-12_dp)
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
  
  call get_lattice_pointGroup(R, pgOps, eps)
  
  ! Normal tests
  call pysave(K, "../tests/rhombohedral/K.in.vasp1")
  call pysave(R, "../tests/rhombohedral/R.in.vasp1")
  call pysave(shift, "../tests/rhombohedral/shift.in.vasp1")
  call pysave(klist, "../tests/rhombohedral/unreduced_klist.in.vasp1")
  call pysave(pgOps, "../tests/rhombohedral/symops.in.vasp1")  
  call symmetryReduceKpointList(K, R, shift,  klist, pgOps, rdKlist, weights, eps)
  call pysave(rdKlist, "../tests/rhombohedral/rhombohedral_kpts.out.vasp1")
  call pysave(weights, "../tests/rhombohedral/rhombohedral_wts.out.vasp1")

  ! VASP
  ! call pysave(K, "../tests/tetragonal/K.in.1")
  ! call pysave(R, "../tests/tetragonal/R.in.1")
  ! call pysave(shift, "../tests/tetragonal/shift.in.1")
  ! call pysave(klist, "../tests/tetragonal/unreduced_klist.in.1")
  ! call pysave(pgOps, "../tests/tetragonal/symops.in.1")
  ! call symmetryReduceKpointList(K, R, shift,  klist, pgOps, rdKlist, weights, &
  !      , eps_=eps)
           
  write(*,'(//"**********")')
  
  do i = 1,size(weights)
     write(*,'(3(1x,f20.3),3x,"w:",i5)') rdKlist(i,:),weights(i)
  end do

  write(*,'(//)')
  write(*,'("Unrd kpts: ",i7)') size(klist,1)
  write(*,'("Rdcd kpts: ",i7)') size(rdKlist,1)
  write(*,'("Rdn ratio: ",f6.3)') size(klist,1)/real(size(weights))

END PROGRAM kpoint_driver
