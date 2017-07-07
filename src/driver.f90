PROGRAM kpoint_driver
  use kpointGeneration
  use num_types
  use vector_matrix_utilities
  use symmetry, only : get_lattice_pointGroup
  use fortpy, only : pysave
  implicit none

  real(dp) :: K(3,3), R(3,3), Hinv(3,3), eps, shift(3)
  integer H(3,3),i
  real(dp), allocatable :: klist(:,:)
  real(dp), pointer     :: pgOps(:,:,:), rdKlist(:,:)
  integer, pointer  :: weights(:)
  eps = 1e-10_dp

  H = transpose(reshape((/ 3, 0, 0, &
                           1, 2, 0, &
                           0, 0, 1/),(/3,3/))) !BCC test 6  
  
  ! R = transpose(reshape((/1.0_dp, 0.0_dp, 0.0_dp, &
  !                         0.0_dp,  1.0_dp, 0.0_dp, &
  !                         0.0_dp,  0.0_dp, 1.0_dp/),(/3,3/)))

  R = transpose(reshape((/-0.4_dp, 0.4_dp, 0.4_dp, &
                          0.4_dp,  -0.4_dp, 0.4_dp, &
                          0.4_dp,  0.4_dp, -0.4_dp/),(/3,3/))) !BCC test 6
  
  call matrix_inverse(real(H,dp),Hinv,eps_=1e-12_dp)

  K = matmul(R,Hinv)

  write(*,'(3("R: ",3(1x,f7.3),/))') (R(i,:),i=1,3)
  write(*,'(3("K: ",3(1x,f7.3),/))') (K(i,:),i=1,3)
  write(*,'(3("Hinv: ",3(1x,f7.3),/))') (Hinv(i,:),i=1,3)
  write(*,'(3("H: ",3(1x,i3),/))') (H(i,:),i=1,3)

  ! write(*,'(3("PP: ",3(1x,f7.3),/))') matmul(K,(/1,0,0/))
  shift = (/ 0.00001_dp, 0.00001_dp, 0.00001_dp /) !BCC test 6
  ! shift  =  (/1.4_dp,1.4_dp,1.4_dp/) !test 7
  write(*,'("new shift: ",3(f6.3,1x))') matmul(K,shift)
  
  call generateFullKpointList(K, R, shift, klist)

  do i = 1,determinant(H)
     write(*,'(3(1x,g11.4))') klist(i,:)
  end do

  call get_lattice_pointGroup(R, pgOps, eps)

  call pysave(K, "../tests/body_centered_cubic/K.in.6")
  call pysave(R, "../tests/body_centered_cubic/R.in.6")
  call pysave(shift, "../tests/body_centered_cubic/shift.in.6")
  call pysave(klist, "../tests/body_centered_cubic/unreduced_klist.in.6")
  call pysave(pgOps, "../tests/body_centered_cubic/symops.in.6")
  call symmetryReduceKpointList(K, R, shift,  klist, pgOps, rdKlist, weights, eps)
  call pysave(rdKlist, "../tests/body_centered_cubic/simple_cubic_kpts.out.6")
  call pysave(weights, "../tests/body_centered_cubic/simple_cubic_wts.out.6")

   write(*,'(//"**********")')

   do i = 1,size(weights)
      write(*,'(3(1x,f6.3),3x,"w:",i5)') rdKlist(i,:),weights(i)
   end do

  write(*,'(//)')
  write(*,'("Unrd kpts: ",i7)') size(klist,1)
  write(*,'("Rdcd kpts: ",i7)') size(rdKlist,1)
  write(*,'("Rdn ratio: ",f6.3)') size(klist,1)/real(size(weights))

END PROGRAM kpoint_driver
