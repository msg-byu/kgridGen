PROGRAM kpoint_driver
  use kpointGeneration
  use num_types
  use vector_matrix_utilities
  use symmetry, only : get_lattice_pointGroup
  implicit none
  
  real(dp) K(3,3), R(3,3), Hinv(3,3), eps
  integer H(3,3),i
  real(dp), allocatable :: klist(:,:)
  real(dp), pointer     :: pgOps(:,:,:), rdKlist(:,:)
  integer, pointer  :: weights(:)
  
  eps = 1e-10_dp
  H = transpose(reshape((/ 2, 0, 0, &
                 1, 3, 0, &
                 0, 0, 1/),(/3,3/)))


  R = transpose(reshape((/-0.5_dp, 0.5_dp, 0.5_dp, &
                 0.5_dp,-0.5_dp, 0.5_dp, &
                 0.5_dp, 0.5_dp,-0.5_dp/),(/3,3/)))
  R = transpose(reshape((/1.0_dp, 0.0_dp, 0.0_dp, &
                 0.0_dp,  1.0_dp, 0.0_dp, &
                 0.0_dp,  0.0_dp, 1.0_dp/),(/3,3/)))

  call matrix_inverse(real(H,dp),Hinv,eps_=1e-12_dp)

  K = matmul(R,Hinv)

  write(*,'(3("R: ",3(1x,f7.3),/))') (R(i,:),i=1,3)
  write(*,'(3("K: ",3(1x,f7.3),/))') (K(i,:),i=1,3)
  write(*,'(3("Hinv: ",3(1x,f7.3),/))') (Hinv(i,:),i=1,3)
  write(*,'(3("H: ",3(1x,i3),/))') (H(i,:),i=1,3)

!  write(*,'(3("PP: ",3(1x,f7.3),/))') matmul(K,(/1,0,0/))
  
  call generateFullKpointList(K, R, klist)

  do i = 1,determinant(H)
     write(*,'(3(1x,g11.4))') klist(i,:)
  end do

  call get_lattice_pointGroup(K, pgOps, eps)
  call symmetryReduceKpointList(K, R, klist, pgOps, rdKlist, weights, eps)


END PROGRAM kpoint_driver

  
