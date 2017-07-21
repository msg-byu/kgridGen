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
  eps = 1e-10_dp
  
  ! Reciprocal lattice vectors
  ! R = transpose(reshape((/ 1.98520863_dp,    0.00000000_dp,   0.00000000_dp, &
  !                          0.0_dp,           1.44640546_dp,   0.00000000_dp, &
  !                          0.0575324872_dp,  0.0_dp,          1.42600347_dp /),(/3,3/)))

  ! R = transpose(reshape((/     /),(/3,3/)))
  R = transpose(reshape((/ 2.222_dp, 0.0_dp, 0.0_dp, &
                           0.0_dp, 2.8222_dp, 0.0_dp, &
                           0.0_dp, 0.0_dp, 2.222o_dp  /),(/3,3/)))
  
  ! HNF Matrix
  H = transpose(reshape((/ 41,  0,  0, &
                            0,   21,  0, &
                            0,   0,   41 /),(/3,3/)))
  
  shift = (/ 0.0_dp, 0.0_dp, 0.0_dp /)
  
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
  ! call pysave(K, "../tests/simple_cubic/K.in.10")
  ! call pysave(R, "../tests/simple_cubic/R.in.10")
  ! call pysave(shift, "../tests/simple_cubic/shift.in.10")
  ! call pysave(klist, "../tests/simple_cubic/unreduced_klist.in.10")
  ! call pysave(pgOps, "../tests/simple_cubic/symops.in.10")  
  ! call symmetryReduceKpointList(K, R, shift,  klist, pgOps, rdKlist, weights, eps)
  ! call pysave(rdKlist, "../tests/simple_cubic/simple_cubic_kpts.out.10")
  ! call pysave(weights, "../tests/simple_cubic/tetragonal_wts.out.10")

  ! VASP
  ! call pysave(K, "../tests/body-centered_tetragonal/K.in.vasp10")
  ! call pysave(R, "../tests/body-centered_tetragonal/R.in.vasp10")
  ! call pysave(shift, "../tests/body-centered_tetragonal/shift.in.vasp10")
  ! call pysave(klist, "../tests/body-centered_tetragonal/unreduced_klist.in.vasp10")
  ! call pysave(pgOps, "../tests/body-centered_tetragonal/symops.in.vasp10")
  call symmetryReduceKpointList(K, R, shift,  klist, pgOps, rdKlist, weights, &
       eps_=eps)
  ! call pysave(rdKlist, "../tests/body-centered_tetragonal/body-centered_tetragonal_kpts.out.vasp10")
  ! call pysave(weights, "../tests/body-centered_tetragonal/body-centered_tetragonal_wts.out.vasp10")
           
  write(*,'(//"**********")')
  
  do i = 1,size(weights)
     write(*,'(3(1x,f20.3),3x,"w:",i5)') rdKlist(i,:),weights(i)
  end do

  write(*,'(//)')
  write(*,'("Unrd kpts: ",i7)') size(klist,1)
  write(*,'("Rdcd kpts: ",i7)') size(rdKlist,1)
  write(*,'("Rdn ratio: ",f6.3)') size(klist,1)/real(size(weights))

END PROGRAM kpoint_driver
