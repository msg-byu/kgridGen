PROGRAM kpoint_driver
  use kpointGeneration
  use num_types
  use vector_matrix_utilities
  use symmetry, only : get_lattice_pointGroup
  ! use fortpy, only : pysave
  implicit none

  real(dp) :: K(3,3), R(3,3), Hinv(3,3), eps, shift(3)
  integer H(3,3),i
  real(dp), allocatable :: klist(:,:)
  real(dp), pointer     :: pgOps(:,:,:), rdKlist(:,:)
  integer, pointer  :: weights(:)
  eps = 1e-10_dp

  H = transpose(reshape((/ 2, 0, 0, &
                           0, 2, 0, &
                           0, 0, 2/),(/3,3/))) !test 11

  ! R = transpose(reshape((/-0.5_dp, 0.5_dp, 0.5_dp, &
  !                0.5_dp,-0.5_dp, 0.5_dp, &
  !                0.5_dp, 0.5_dp,-0.5_dp/),(/3,3/)))
  
  R = transpose(reshape((/1.0_dp, 0.0_dp, 0.0_dp, &
                          0.0_dp,  1.0_dp, 0.0_dp, &
                          0.0_dp,  0.0_dp, 1.0_dp/),(/3,3/)))

  call matrix_inverse(real(H,dp),Hinv,eps_=1e-12_dp)

  K = matmul(R,Hinv)

  write(*,'(3("R: ",3(1x,f7.3),/))') (R(i,:),i=1,3)
  write(*,'(3("K: ",3(1x,f7.3),/))') (K(i,:),i=1,3)
  write(*,'(3("Hinv: ",3(1x,f7.3),/))') (Hinv(i,:),i=1,3)
  write(*,'(3("H: ",3(1x,i3),/))') (H(i,:),i=1,3)

  ! write(*,'(3("PP: ",3(1x,f7.3),/))') matmul(K,(/1,0,0/))

  shift  =  (/1.0_dp/3.0_dp, 1.0_dp/3.0_dp, 1.0_dp/3.0_dp/)
  call generateFullKpointList(K, R, shift, klist)

  ! do i = 1,determinant(H)
  !    write(*,'(3(1x,g11.4))') klist(i,:)
  ! end do

  call get_lattice_pointGroup(K, pgOps, eps)

  ! call pysave(K, "../tests/simple_cubic/K.in.2")
  ! call pysave(R, "../tests/simple_cubic/R.in.2")
  ! call pysave(shift, "../tests/simple_cubic/shift.in.2")
  ! call pysave(klist, "../tests/simple_cubic/unreduced_klist.in.2")
  ! call pysave(pgOps, "../tests/simple_cubic/symops.in.2")
  call symmetryReduceKpointList(K, R, shift,  klist, pgOps, rdKlist, weights, eps)
  ! call pysave(rdKlist, "../tests/simple_cubic/simple_cubic_kpts.out.2")
  ! call pysave(weights, "../tests/simple_cubic/simple_cubic_wts.out.2")

   write(*,'(//"**********")')

   do i = 1,size(weights)
      write(*,'(3(1x,f6.3),3x,"w:",i5)') rdKlist(i,:),weights(i)
   end do

  write(*,'(//)')
  write(*,'("Unrd kpts: ",i7)') size(klist,1)
  write(*,'("Rdcd kpts: ",i7)') size(rdKlist,1)
  write(*,'("Rdn ratio: ",f6.3)') size(klist,1)/real(size(weights))

END PROGRAM kpoint_driver
