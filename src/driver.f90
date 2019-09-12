PROGRAM kpoint_driver
  USE kpointGeneration
  USE num_types
  USE vector_matrix_utilities
  USE symmetry, ONLY : get_lattice_pointGroup
  USE rational_mathematics, ONLY: HermiteNormalForm
  IMPLICIT NONE
  REAL(dp)              :: K(3,3), R(3,3), Hinv(3,3), reps, aeps, shift(3)
  REAL(dp)              :: H(3,3)
  REAL(dp), POINTER     :: klist(:,:)
  REAL(dp), POINTER     :: pgOps(:,:,:), rdKlist(:,:)
  INTEGER, POINTER      :: weights(:)
  INTEGER i
  ! integer H(3,3), i

  ! shift = (/ 0.0_dp, 0.0_dp, 0.0_dp /)
  ! R = transpose(reshape((/ -0.0_dp, 1.0_dp, 1.0_dp, &
  !                         1.0_dp, 0.0_dp, 1.0_dp, &
  !                         1.0_dp, 1.0_dp, 0.0_dp /),(/3,3/)))
  ! H = transpose(reshape((/ 2, 0, 0, &
  !                         0, 2, 0, &
  !                         0, 0, 2 /),(/3,3/)))

  ! Finite precision tolerance (same as default value)
  ! eps = 1e-10_dp
  ! H = real(H,dp)
  ! call matrix_inverse(real(H,dp), Hinv, eps_=eps)

  ! Columns of K are the grid generating vectors.
  ! K = matmul(R,Hinv)


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

  ! write(*,'(3("R: ",3(1x,f7.3),/))') (R(i,:),i=1,3)
  ! write(*,'(3("H: ",3(1x,i3),/))') (H(i,:),i=1,3)
  ! write(*,'(3("Hinv: ",3(1x,f7.3),/))') (Hinv(i,:),i=1,3)
  ! write(*,'(3("K: ",3(1x,f11.7),/))') (K(i,:),i=1,3)
  ! write(*,'("shift: ",3(f6.3,1x))') shift
  ! write(*,'("cart shift: ",3(f6.3,1x))') matmul(K,shift)

  ! write(*,'(3("PP: ",3(1x,f7.3),/))') matmul(K,(/1,0,0/))
  ! call generateFullKpointList(K, R, shift, klist, eps)
  ! do i = 1,determinant(H)
  !    write(*,'(3(1x,g11.4))') klist(i,:)
  ! end do
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

  ! write(*,'(//)')
  ! write(*,'("Unrd kpts: ",i7)') size(klist,1)
  ! write(*,'("Rdcd kpts: ",i7)') size(rdKlist,1)
  ! write(*,'("Rdn ratio: ",3x,f4.1)') size(klist,1)/real(size(weights))


  ! Map into first Brillouin zone tests
  reps = 1e-8_dp
  aeps = 1e-10_dp

  shift = (/ 0.0_dp, -1.5_dp, 0.5_dp /)
  R = TRANSPOSE(RESHAPE((/ 0.0_dp, 0.86172861889_dp, 0.86172861889_dp, &
       0.86172861889_dp, 0.0_dp, 0.86172861889_dp, &
       0.86172861889_dp, 0.86172861889_dp, 0.0_dp /),(/3,3/)))
  H = TRANSPOSE(RESHAPE((/ 2, 0, 0, &
       2, 3, 0, &
       3, 1, 4 /),(/3,3/)))

  CALL matrix_inverse(REAL(H,dp), Hinv, eps_=aeps)
  ! Columns of K are the grid generating vectors.
  K = MATMUL(R,Hinv)
  WRITE(*,'(3("K: ",3(1x,f11.7),/))') (K(i,:),i=1,3)
  CALL get_lattice_pointGroup(R, pgOps, eps_=aeps)
  CALL generateFullKpointList(K, R, shift, klist, reps_=reps,aeps_=aeps)
  WRITE(*,'("shift: ",3(f6.3,1x))') shift
  CALL symmetryReduceKpointList(K, R, shift,  klist, pgOps, rdKlist, weights, &
       reps_=reps,aeps_=aeps)


  WRITE(*,'(//"**********")')
  DO i = 1,SIZE(rdKlist,1)
     WRITE(*,'(3(1x,f9.3),3x,"w:",i5)') rdKlist(i,:),weights(i)
  END DO

  WRITE(*,'(//)')
  WRITE(*,'("Unrd kpts: ",i7)') SIZE(klist,1)
  WRITE(*,'("Rdcd kpts: ",i7)') SIZE(rdKlist,1)
  WRITE(*,'("Rdn ratio: ",3x,f4.1)') SIZE(klist,1)/REAL(SIZE(weights))

  CALL mapKptsIntoBZ(R, rdKlist, reps_=reps, aeps_=aeps)

  WRITE(*,'(//"**********")')
  DO i = 1,SIZE(rdKlist,1)
     WRITE(*,'(3(1x,f9.3),3x,"w:",i5)') rdKlist(i,:),weights(i)
  END DO
  WRITE(*,'(//"**********")')


END PROGRAM kpoint_driver
