MODULE kpointGeneration
  use num_types
  use numerical_utilities, only: equal
  use vector_matrix_utilities, only: matrix_inverse, norm, &
       cross_product, volume, determinant
  use rational_mathematics, only: HermiteNormalForm, SmithNormalForm
  use symmetry
  implicit none

  private
  public generateFullKpointList, symmetryReduceKpointList
CONTAINS

  !!<summary>Takes two lattices, a generating lattice (K) and the
  !! reciprocal lattice (R), and generates all the points of K that
  !! are inside one unit cell of R. </summary>
  !!<parameter name="K" regular="true">3x3 matrix. Columns are the
  !! generating vectors of the k-grid</parameter>
  !!<parameter name="R" regular="true">3x3 matrix. Columns are the
  !! reciprocal lattice vectors</parameter>
  !!<parameter name="shift" regular="true">Fractional shift of the
  !! k-grid. Given as fractions of the generating k-lattice vectors.
  !!</parameter> 
  !!<parameter name="KpList" regular="true">Full list of (unreduced)
  !! k-points</parameter>  
  subroutine generateFullKpointList(K, R, shift, KpList, eps_)
    real(dp), intent(in) :: K(3,3)
    real(dp), intent(in) :: R(3,3)
    real(dp), intent(in) :: shift(3)
    real(dp), allocatable:: KpList(:,:)
    real(dp), intent(in), optional:: eps_   

    ! Inverse of the kgrid matrix, Inverse of reciprocal cell, the
    ! shift in Cartesian coordinates
    real(dp) :: Kinv(3,3), Rinv(3,3), carShift(3)
    ! Finite precision parameter
    real(dp) :: eps 
    ! Integer matrices for HNF conversion
    integer  :: S(3,3), H(3,3), B(3,3)
    ! Number of kpoints (i.e., index of kgrid lattice, a.k.a. volume
    ! factor)
    integer  :: n
    ! Diagonal entries of the HNF matrix
    integer  :: a, c, f
    ! loop counters of kpoint generating vectors
    integer  :: iK, jK, kK
    ! loop over kpoints
    integer  :: iKP
    ! index (ordinal number) of kpoint
    integer  :: idx
    ! flag for catching errors
    logical  :: err
    
    if(.not. present(eps_)) then
       eps = 1e-10_dp
    else
       eps =  eps_
    endif

    ! Check for valid inputs
    if (determinant(K) > determinant(R)+eps) then
       write(*,*) "ERROR (generateFullKpointList in &
            &generateKpoints.f90):"
       stop "The kpoint generating lattice vectors have a unit cell&
            & larger than the reciprocal lattice."
    endif
    call matrix_inverse(K,Kinv,err)
    if (err) then
       write(*,*) "ERROR: (generateFullKpointList in &
            &generateKpoints.f90):"
       write(*,*) "The kgrid generating vectors are linearly &
            &dependent."
       stop
    end if
    if (any(matmul(Kinv,R) -  nint(matmul(Kinv,R)) > eps)) then
       write(*,*) "ERROR: (generateFullKpointList in &
            &generateKpoints.f90):"
       stop "The point generating vectors and the reciprocal lattice &
            &are incommensurate."
    endif
    call matrix_inverse(R,Rinv,err)
        if (err) then
       write(*,*) "ERROR: (generateFullKpointList in &
            &generateKpoints.f90):"
       write(*,*) "The reciprocal lattice vectors are linearly &
            &dependent."
       stop
    end if
    if (any(matmul(Kinv,R) -  nint(matmul(Kinv,R)) > eps)) then
       write(*,*) "ERROR: (generateFullKpointList in &
            &generateKpoints.f90):"
       stop "The point generating vectors and the reciprocal lattice &
            &are incommensurate."
    endif

    ! Put the Shift in Cartesian coordinates.
    carShift = matmul(K,shift)
    ! Bring the shift back into the first unit cell.
    call bring_into_cell(carShift, Rinv, R, eps)
    
    ! Integer transformation matrix that takes K to R, not necessarily
    ! HNF at the outset
    S = nint(matmul(Kinv,R))
    ! Find the HNF of S, store it in H (B is the transformation
    ! matrix)
    call HermiteNormalForm(S,H,B)
    ! This the volume ratio between R and K, and therefore the number
    ! of unreduced kpoints
    n = determinant(H) 
    a = H(1,1); c = H(2,2); f = H(3,3)
    allocate(KpList(n,3))
    ! Generate a list of kpoints (not necessarily all in one unit
    ! cell) but which are unique under lattice translations.
    do iK = 0, a-1
       do jK = 0, c-1
          do kK = 0, f-1
             ! ordinal counter for k-points;
             ! converted to base-10 from mixed-radix number of HNF
             ! diagonal entries
             idx = f*c*iK + f*jK + kK + 1
             ! compute Cartesian coordinates of the k-point
             KpList(idx,:) = matmul(K,(/iK, jK, kK/)) + carShift
          enddo
       enddo
    enddo

    ! do idx = 1,n
    !    write(*,'("kp: ",3(1x,f8.3))') KpList(idx,:)
    ! end do
    
    ! Bring each k-point into the first unit cell
    do iKP = 1,n
       call bring_into_cell(KpList(iKP,:),Rinv,R,eps)
    enddo


    ! do idx = 1,n
    !    write(*,'("kp (bic): ",3(1x,f8.3))') KpList(idx,:)
    ! end do

  END subroutine generateFullKpointList

  !!<summary>Reduce k-points to irreducible set. Takes a list of
  !! translationally distinct, but not rotationally distinct, kpoints
  !! and reduces them by the given point group. </summary>
  !!<parameter regular="true" name="K">Matrix of grid generating
  !! vectors.</parameter>
  !!<parameter regular="true" name="R">Matrix of reciprocal lattice
  !! vectors.</parameter>
  !!<parameter regular="true" name="shift">Offset (shift) of the
  !! k-grid (fractions of k-grid vectors).</parameter>
  !!<parameter regular="true" name="UnreducedKpList">List of
  !! unreduced kpoints.</parameter>
  !!<parameter regular="true" name="SymOps">Integer rotation matrices
  !! (coordinates of the lattice basis).</parameter>
  !!<parameter regular="true" name="ReducedList">List of symmetry-
  !! reduced kpoints.</parameter>
  !!<parameter regular="true" name="weights"> "Weights" of kpoints
  !! (length of each orbit).</parameter>
  !!<parameter regular="true" name="eps_">Finite precision parameter
  !! (optional)</parameter>
  subroutine symmetryReduceKpointList(K, R, shift, UnreducedKpList, &
       SymOps, ReducedList, weights, eps_)
    ! basis vectors of the k-point grid, basis vectors of recip.
    ! lattice
    real(dp), intent(in) :: K(3,3), R(3,3)
    ! Offset of the k-point grid
    real(dp), intent(in) :: shift(3) 
    ! kpoint number (row); 3 components (columns)
    real(dp), intent(in) :: UnreducedKpList(:,:)
    ! Last slot is op # index, first two are 3x3
    real(dp), intent(in) :: SymOps(:,:,:)
    real(dp), pointer    :: ReducedList(:,:)
    integer, pointer     :: weights(:) 
    real(dp), optional   :: eps_

    integer :: iOp, nOps, iRdKpt, nRdKpt, iUnRdKpt, nUR, cOrbit
    integer :: idx, i, sum, hashTable(size(UnreducedKpList,1))
    ! Inverse of kgrid matrix, Inverse of reciprocal lattice
    real(dp):: InvK(3,3), InvR(3,3)
    ! unrotated kpoint, rotated kpoint, kpoint in gspace coords
    real(dp):: urKpt(3), roKpt(3), gpt(3) 
    ! coords (unrotated, rotated), in g-space
    ! Integer transformation that takes K to R
    integer :: N(3,3) 
    ! HNF, HNF transform matrix, SNF transform matrices, SNF, diag(SNF)
    integer :: H(3,3), B(3,3), L(3,3), Ri(3,3), S(3,3), D(3) 
    real(dp):: eps
    logical :: err
    
    if(.not. present(eps_)) then
       eps = 1e-10_dp
    else
       eps =  eps_
    endif

    call matrix_inverse(K,InvK,err,eps)
    if (err) then
       write(*,*) "ERROR: (symmetryReduceKpointList in &
            &generateKpoints.f90)"
       stop "The k-grid  vectors that were passed in are not linearly&
            &dependent."
    END if

    call matrix_inverse(R,InvR,err,eps)
    if (err) then
       write(*,*) "ERROR: (symmetryReduceKpointList in &
            &generateKpoints.f90)"
       stop "The reciprocal lattice vectors that were passed in are &
            &not linearly dependent."
    END if
    
    !! Check for valid inputs
    ! Make sure kgrid is commensurate with reciprocal cell
    ! if (any(matmul(InvK,R) - nint(matmul(InvK,R)) > eps)) then
    if (.not. equal(matmul(InvK, R), nint(matmul(InvK, R)), eps)) then
       write(*,*) "ERROR: (symmetryReduceKpointList in &
            &generateKpoints.f90):"
       write(*,*) "The point generating vectors and the reciprocal &
            &lattice are incommensurate."
       write(*,*) "The reciprocal cell should be an (matrix) integer &
            &multiple of the generating lattice vectors."
       stop
    endif
    ! Make sure that kgrid is no bigger than reciprocal cell
    if (determinant(K) > determinant(R)+eps) then
       write(*,*) "ERROR (symmetryReduceKpointList in &
            &generateKpoints.f90):"
       write(*,*) "The kpoint generating lattice vectors define a &
            &unit cell larger than the reciprocal lattice."
       write(*,*) "This doesn't make sense. There should be at least &
            &1 kpoint per reciprocal cell."
       stop
    endif

    ! Integer transformation matrix that takes K to R, not necessarily
    ! HNF at the outset
    N = nint(matmul(InvK,R))
    ! Find the HNF of N, store it in H (B is the transformation
    ! matrix)
    call HermiteNormalForm(N,H,B)
    ! Left side of transform will be used later, right side (Ri)
    ! will not be.
    call SmithNormalForm(H,L,S,Ri)
    write(*,'(3(1x,i3))') (H(i,:),i=1,3)
    write(*,*) 
    write(*,'(3(1x,i3))') (S(i,:),i=1,3)
    write(*,'("InvK")')
    write(*,'(3(1x,f7.3))') (InvK(i,:),i=1,3)
    
    D = (/S(1,1),S(2,2),S(3,3)/)
    
    ! Count how many unique orbits there are; i.e., the number of
    ! unique kpoints
    cOrbit = 0
    ! Number of symmetry operators in the point group
    nOps = size(SymOps,3)
    ! Number of unreduced kpoints
    nUR = size(UnreducedKpList,1)
    ! Use a hash table to keep track of symmetrically-equivalent
    ! k-points
    hashTable = 0
    ! Loop over each kpoint and mark off its symmetry brothers
    do iUnRdKpt = 1,nUR
       ! unrotated kpoint (shorter name for clarity's sake)
       ! urKpt = UnreducedKpList(iUnRdKpt,:) 
       urKpt = UnreducedKpList(iUnRdKpt,:)
       write(*,'(/,"kpt#: ",i4)') iUnRdKpt
       idx = findKptIndex(urKpt,K,InvK,L,D, shift)
       write(*,'("urKpt: ",3(f6.3,1x),"idx:",i3)') urKpt, idx
       ! This kpoint is already on an orbit, skip it
       if (hashTable(idx)/=0) cycle
       cOrbit = cOrbit + 1
       hashTable(idx) = cOrbit
       ! Loop over each symmetry operator, applying it to each
       ! k-point
       do iOp = 1,nOps
          ! Rotate the kpoint
          roKpt = matmul(SymOps(:,:,iOp),urKpt)
          write(*,'(/,"Operator:")') 
          write(*,'(3(1x,f7.3))') (SymOps(i,:,iOp),i=1,3)
          write(*,'("roKpt: ",3(f6.3,1x),i3)') roKpt
          ! Make sure the rotated kpoint is still part of the k-grid.
          ! If not, cycle
          call bring_into_cell(roKpt, InvR, R, eps)
          write(*,'("icKpt: ",3(f6.3,1x),i3)') roKpt
          ! Unshift the k-point temporarily to check if it still is
          ! a member of the lattice
          roKpt = roKpt - matmul(K,shift)
          ! write(*,'(3(1x,i3))') nint(matmul(invK,roKpt))
          ! write(*,'(3(1x,f7.3))') (matmul(invK,roKpt) - nint(matmul(invK,roKpt)))
          if (.not. equal(matmul(invK,roKpt), &
               nint(matmul(invK,roKpt)), eps)) then
             cycle
          endif
          ! Reshift the kpoint before finding its indx
          roKpt = roKpt + matmul(K,shift)
          write(*,'("shKpt: ",3(f6.3,1x),i3)') roKpt
          idx = findKptIndex(roKpt, K, InvK, L, D, shift)
          write(*,'("Op#",1x,i2,5x,"rkpt: &
               &",3(f6.3,1x),5x,"idx: ",i4)') iOp,roKpt,idx
          hashTable(idx) = cOrbit
       enddo
    enddo
    write(*,*) "Hash table:"
    do i = 1, nUr
       write(*,'("kpt#:",i3,3x,"index:",i3)') i, hashTable(i)
    enddo
    ! Now that we have the hash table populated, make a list of the
    ! irreducible kpoints and their corresponding weights.
    sum = 0
    allocate(ReducedList(cOrbit,3),weights(cOrbit))
    do i = 1, cOrbit
       weights(i) = count(hashTable==i)
       sum = sum + weights(i)
       ! A hack to get the location of the first match
       idx = minloc(hashTable, dim=1, mask=hashTable==i)
       ! print*,"idx",idx
       ReducedList(i,:) = UnreducedKpList(idx,:)
    enddo    
    ! Fail safe checks
    if (any(hashTable==0)) then
       write(*,*) "ERROR: (symmetryReduceKpointList in &
            &generateKpoints.f90)"
       write(*,*) "At least one k-point in the unreduced list was &
            &not included in one of the orbits of the &
            &symmetry group."
       stop
    endif
    if (size(UnreducedKpList,1)/real(size(weights)) > size(SymOps,3)) then
       write(*,*) "ERROR: (symmetryReduceKpointList in &
            &generateKpoints.f90)"
       write(*,*) "The ratio of unreduced to reduced k-points is larger than &
            &the size of the point group."
       stop
    endif
    if (.not. sum == nUr) then
       write(*,*) "ERROR: (symmetryReduceKpointList in &
            &generateKpoints.f90)"
       write(*,*) "The sum of the weights of the irreducible k-points&
            & is not equal to the number of unreduced k-points."
    endif
    
  CONTAINS
    ! This function takes a k-point in Cartesian coordinates and
    ! "hashes" it into a single number, corresponding to its
    ! place in the k-point list.
    function findKptIndex(kpt, K, InvK, L, D, shift)
      ! The k-point, k-grid generating vecs, inverse of the k-grid
      ! generating vecs
      real(dp), intent(in) :: kpt(3), K(3,3), InvK(3,3), shift(3)
      ! Left transform for SNF conversion, diagonal of SNF
      integer,  intent(in) :: L(3,3), D(3)
      ! index of the kpoint (base 10, 1..n)
      integer              :: findKptIndex
      ! Fractional shift of the k-grid given as fractions of the
      ! generating k-lattice vectors, k-point in lattice coordinates,
      ! temporary k-point without shift
      real(dp) :: gpt(3), tmpKpt(3)
    
      ! Unshift the k-point temporarily to check if it still is a
      ! member of the lattice
      tmpKpt = kpt - matmul(K,shift)
      ! write(*,'(3(1x,i3))') nint(matmul(invK,roKpt))
      ! write(*,'(3(1x,f7.3))') (matmul(invK,roKpt) - &
      !      nint(matmul(invK,roKpt)))
      ! if (.not. equal(matmul(invK,roKpt), &
      !      nint(matmul(invK,roKpt)), eps)) then
      ! if (any((matmul(invK,roKpt) - &
      !      nint(matmul(invK,roKpt))) > eps)) then
      if (.not. equal(matmul(invK,tmpKpt), nint(matmul(invK,tmpKpt)),&
           eps)) then
         write(*,*) "ERROR: (symmetryReduceKpointList in &
              &kpointGeneration)"
         write(*,*) "A rotated k-point is not a lattice point of the &
              &generating vectors."
         stop
      endif
            
      ! k-point in lattice coordinates
      gpt = matmul(InvK,kpt)
      ! Convert the kpoint to group coordinates and bring into first
      ! unit cell
      gpt = modulo(matmul(L,nint(gpt)),D)
      ! Convert from group coordinates (3-vector) to single base-10
      ! number between 1 and nUr
      ! Hash of the kpt
      findKptIndex = gpt(1)*D(2)*D(3) + gpt(2)*D(3) + gpt(3) + 1  
    END function findKptIndex
    
  END subroutine symmetryReduceKpointList
  
END MODULE kpointGeneration

