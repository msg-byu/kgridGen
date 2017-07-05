MODULE kpointGeneration
  use num_types
  use numerical_utilities, only: equal
  use vector_matrix_utilities, only: matrix_inverse, norm, cross_product, volume, &
       determinant
  use rational_mathematics, only: HermiteNormalForm, SmithNormalForm
  use symmetry
  implicit none

  private
  public generateFullKpointList, symmetryReduceKpointList
CONTAINS

  !!<summary>Takes two lattices, a generating lattice (K) and the reciprocal lattice (R),
  !! and generates all the points of K that are inside one unit cell of R. </summary>
  !!<parameter name="K" regular="true">3x3 matrix. Columns are the generating vectors of
  !! the k-grid</parameter>
  !!<parameter name="R" regular="true">3x3 matrix. Columns are the reciprocal lattice
  !! vectors</parameter>
  !!<parameter name="kLVshift" regular="true">Fractional shift of the k-grid. Given as
  !! fractions of the generating k-lattice vectors. </parameter> 
  !!<parameter name="KpList" regular="true">Full list of (unreduced) k-points</parameter>  
  subroutine generateFullKpointList(K, R, kLVshift, KpList, eps_)
    real(dp), intent(in) :: K(3,3)
    real(dp), intent(in) :: R(3,3)
    real(dp), intent(in) :: kLVshift(3)
    real(dp), allocatable:: KpList(:,:)
    real(dp), intent(in), optional:: eps_   

    real(dp) :: Kinv(3,3) ! Inverse of the k-grid cell 
    real(dp) :: Rinv(3,3) ! Inverse of reciprocal cell
    real(dp) :: eps ! Finite precision parameter
    real(dp) :: cartShift(3), bicCartShift(3) ! k-grid shifts in Cartesian coordinates
    integer  :: S(3,3), H(3,3), B(3,3) ! Integer matrices for HNF conversion
    integer  :: n ! Number of k-points (i.e., index of kgrid lattice, a.k.a. volume factor)
    integer  :: a, c, f ! Diagonal entries of the HNF matrix
    integer  :: iK, jK, kK ! loop counters of k-point generating vectors
    integer  :: iKP ! loop over k-points
    integer  :: idx ! index (ordinal number) of k-point
    logical  :: err ! flag for catching errors
    
    if(.not. present(eps_)) then
       eps = 1e-10_dp
    else
       eps =  eps_
    endif

    ! Check for valid inputs
    if (determinant(K) > determinant(R)+eps) then
       write(*,*) "ERROR (generateFullKpointList in generateKpoints.f90):"
       write(*,*) "The k-point generating lattice vectors have a unit cell &
            &larger than the reciprocal lattice."
       stop
    endif
    call matrix_inverse(K,Kinv,err)
    if (err) then
       write(*,*) "ERROR: (generateFullKpointList in generateKpoints.f90):"
       write(*,*) "The kgrid generating vectors are linearly dependent."
       stop
    endif
    if (any(matmul(Kinv,R) -  nint(matmul(Kinv,R)) > eps)) then
       write(*,*) "ERROR: (generateFullKpointList in generateKpoints.f90):"
       stop "The point generating vectors and the reciprocal lattice are incommensurate."
    endif
    call matrix_inverse(R,Rinv,err)
        if (err) then
       write(*,*) "ERROR: (generateFullKpointList in generateKpoints.f90):"
       write(*,*) "The reciprocal lattice vectors are linearly dependent."
       stop
    endif
    if (any(matmul(Kinv,R) -  nint(matmul(Kinv,R)) > eps)) then
       write(*,*) "ERROR: (generateFullKpointList in generateKpoints.f90):"
       stop "The point generating vectors and the reciprocal lattice are incommensurate."
    endif
    cartShift = matmul(K,kLVshift)
    bicCartShift = cartShift ! The k-grid shift, but bring into cell
    call bring_into_cell(bicCartShift, K, Kinv, eps)
    if (.not. equal(cartShift, bicCartShift, eps)) then
       write(*,*) "WARNING: (generateFullKpointList in generateKpoints.f90)"
       write(*,*) "The given shift was outside the first k-grid cell. By translation"
       write(*,*) "symmetry, this is always equivalent to a shift inside the cell."
    endif
    
    ! Integer transformation matrix that takes K to R, not necessarily HNF at the outset
    S = nint(matmul(Kinv,R))
    ! Find the HNF of S, store it in H (B is the transformation matrix)
    call HermiteNormalForm(S,H,B)
    ! This the volume ratio between R and K, i.e., the number of unreduced k-points
    n = determinant(H) 
    a = H(1,1); c = H(2,2); f = H(3,3)
    allocate(KpList(n,3))
    ! Generate a list of k-points (not necessarily all in one unit cell but)
    ! which are unique under lattice translations
    do iK = 0, a-1
       do jK = 0, c-1
          do kK = 0, f-1
             !ordinal counter for k-points;
             !converted to base-10 from mixed-radix number of HNF diagonal entries
             idx = f*c*iK + f*jK + kK + 1
             ! compute Cartesian coordinates of the k-point
             KpList(idx,:) = matmul(K,(/iK, jK, kK/)) + matmul(K,kLVshift) 
          enddo
       enddo
    enddo

!!do idx = 1,n
!!   write(*,'("kp: ",3(1x,f8.3))') KpList(idx,:)
!!end do
    
    ! Bring each k-point into the first unit cell
    do iKP = 1,n
       call bring_into_cell(KpList(iKP,:),Rinv,R,eps)
    enddo

!!        do idx = 1,n
!!       write(*,'("kp (bic): ",3(1x,f8.3))') KpList(idx,:)
!!    end do

  END subroutine generateFullKpointList

  !!<summary>Reduce k-points to irreducible set. Takes a list of translationally distinct,
  !! but not rotationally distinct, k-points and reduces them by the given point group.
  !! </summary>
  !!<parameter regular="true" name="K">Matrix of grid generating vectors.</parameter>
  !!<parameter regular="true" name="R">Matrix of reciprocal lattice vectors.</parameter>
  !!<parameter regular="true" name="kLVshift">Offset (shift) of the k-grid
  !! (fractions of k-grid vectors).</parameter>
  !!<parameter regular="true" name="UnreducedKpList">Unreduced k-point list.</parameter>
  !!<parameter regular="true" name="SymOps">Integer rotation matrices (coordinates of
  !! the lattice basis).</parameter>
  !!<parameter regular="true" name="ReducedList">List of symmetry-reduced k-points.
  !!</parameter>
  !!<parameter regular="true" name="weights"> "Weights" of k-points (length of each orbit).
  !!</parameter>
  !!<parameter regular="true" name="eps_">Finite precision parameter (optional)</parameter>
  subroutine symmetryReduceKpointList(K, R, kLVshift, UnreducedKpList, SymOps, ReducedList, &
       weights, eps_)
    real(dp), intent(in) :: K(3,3), R(3,3) 
    real(dp), intent(in) :: kLVshift(3) 
    real(dp), intent(in) :: UnreducedKpList(:,:) 
    real(dp), intent(in) :: SymOps(:,:,:) ! Last slot is op # index, first two are 3x3
    real(dp), pointer    :: ReducedList(:,:)
    integer, pointer     :: weights(:) 
    real(dp), optional   :: eps_

    integer :: iOp, nOps, iRdKpt, nRdKpt, iUnRdKpt, nUR, cOrbit, idx, i, sum
    integer :: hashTable(size(UnreducedKpList,1))
    integer :: iFirst(size(UnreducedKpList,1)), iWt(size(UnreducedKpList,1))
    real(dp):: InvK(3,3) ! Inverse of kgrid matrix
    real(dp):: InvR(3,3) ! Inverse of reciprocal lattice
    real(dp):: urKpt(3), roKpt(3) ! unrotated k-point, rotated k-point
    real(dp):: shift(3) ! Shift of k-grid lattice in Cartesian coordinates
    real(dp):: gpt(3) ! k-point in gspace coords  
    integer :: N(3,3) ! Integer transformation that takes K to R
    ! HNF, SNF transform matrices, SNF, diag(SNF)
    integer :: H(3,3), L(3,3), Ri(3,3), S(3,3), D(3) 
    real(dp):: eps
    logical :: err
    
    if(.not. present(eps_)) then
       eps = 1e-10_dp
    else
       eps =  eps_
    endif

    call matrix_inverse(K,InvK,err,eps)
    if (err) then
       write(*,*) "ERROR: (symmetryReduceKpointList in generateKpoints.f90)"
       stop "The k-grid  vectors that were passed in are not linearly dependent."
    END if

    call matrix_inverse(R,InvR,err,eps)
    if (err) then
       write(*,*) "ERROR: (symmetryReduceKpointList in generateKpoints.f90)"
       stop "The reciprocal lattice vectors are not linearly dependent."
    END if

    !! Check for valid inputs
    ! Make sure kgrid is commensurate with reciprocal cell
    if (.not. equal(matmul(InvK,R), nint(matmul(InvK,R)), eps)) then
       write(*,*) "ERROR: (symmetryReduceKpointList in generateKpoints.f90):"
       write(*,*) "The point generating vectors and the reciprocal lattice are &
            &incommensurate."
       write(*,*) "The reciprocal cell should be an (matrix) integer multiple of &
            & the generating lattice vectors."
       stop
    endif
    ! Make sure that kgrid is no bigger than reciprocal cell
    if (abs(determinant(K)) > abs(determinant(R))+eps) then
       write(*,*) "ERROR (symmetryReduceKpointList in generateKpoints.f90):"
       write(*,*) "The k-point generating lattice vectors define a unit cell larger&
            & than the reciprocal lattice."
       write(*,*) "This doesn't make sense. There should be at least 1 k-point per &
            &reciprocal cell."
       stop
    endif

   shift = matmul(K,kLVshift)
   ! Integer transformation matrix that takes K to R, not necessarily HNF at the outset
    N = nint(matmul(InvK,R))
    ! Find the HNF of N, store it in H (B is the transformation matrix)
    call HermiteNormalForm(N,H,L)
    ! Left side of transform will be used later, right side (Ri) will not be.
    call SmithNormalForm(H,L,S,Ri)

    D = (/S(1,1),S(2,2),S(3,3)/)
    
    cOrbit = 0 ! Count number of unique orbits. I.e., the number of unique k-points
    nOps = size(SymOps,3) ! Number of symmetry operators in the point group
    nUR = size(UnreducedKpList,1) ! Number of unreduced k-points
    hashTable = 0 ! Use a hash table to keep track of symmetrically-equivalent k-points
    
    call matrix_inverse(K,InvK,err,eps)
    if (err) then
       write(*,*) "ERROR: (symmetryReduceKpointList in generateKpoints.f90)"
       stop "The k-grid  vectors that were passed in are not linearly dependent."
    END if

    call matrix_inverse(R,InvR,err,eps)
    if (err) then
       write(*,*) "ERROR: (symmetryReduceKpointList in generateKpoints.f90)"
       stop "The reciprocal lattice vectors that were passed in are not linearly dependent."
    END if

    !! Check for valid inputs
    ! Make sure kgrid is commensurate with reciprocal cell
    if (any(matmul(InvK,R) - nint(matmul(InvK,R)) > eps)) then
       write(*,*) "ERROR: (symmetryReduceKpointList in generateKpoints.f90):"
       write(*,*) "The point generating vectors and the reciprocal lattice are incommensurate."
       write(*,*) "The reciprocal cell should be an (matrix) integer multiple of the generating lattice vectors."
       stop
    endif
    ! Make sure that kgrid is no bigger than reciprocal cell
    if (determinant(K) > determinant(R)+eps) then
       write(*,*) "ERROR (symmetryReduceKpointList in generateKpoints.f90):"
       write(*,*) "The kpoint generating lattice vectors define a unit cell larger than the reciprocal lattice."
       write(*,*) "This doesn't make sense. There should be at least 1 kpoint per reciprocal cell."
       stop
    endif

   ! Integer transformation matrix that takes K to R, not necessarily HNF at the outset
    N = nint(matmul(InvK,R))
    ! Find the HNF of N, store it in H (B is the transformation matrix)
    call HermiteNormalForm(N,H,L)
    ! Left side of transform will be used later, right side (Ri) will not be.
    call SmithNormalForm(H,L,S,Ri)
    write(*,'(3(1x,i3))') (H(i,:),i=1,3)
    write(*,*) 
    write(*,'(3(1x,i3))') (S(i,:),i=1,3)
    write(*,'("InvK")')
    write(*,'(3(1x,f7.3))') (InvK(i,:),i=1,3)
    
    D = (/S(1,1),S(2,2),S(3,3)/)
    
    cOrbit = 0 ! Count how many unique orbits there are. I.e., the number of unique kpoints
    nOps = size(SymOps,3) ! Number of symmetry operators in the point group
    nUR = size(UnreducedKpList,1) ! Number of unreduced kpoints
    hashTable = 0 ! Use a hash table to keep track of symmetrically-equivalent kpoints
    iFirst = 0 ! store the index of the first kpoint in an orbit
    iWt = 0 ! count the number of members of each orbit
    do iUnRdKpt = 1,nUR ! Loop over each k-point and mark off its symmetry brothers
       urKpt = UnreducedKpList(iUnRdKpt,:) ! unrotated k-point (shorter name for clarity)
      write(*,'(/,"kpt#: ",i4)') iUnRdKpt
     idx = findKptIndex(urKpt-shift, InvK, L, D)
                write(*,'("urKpt: ",3(f6.3,1x),"idx:",i3)') urKpt, idx
    
     if (hashTable(idx)/=0) cycle ! This k-point is already on an orbit, skip it
     cOrbit = cOrbit + 1
     hashTable(idx) = cOrbit
     iFirst(cOrbit) = idx
     iWt(cOrbit) = 1
     do iOp = 1,nOps ! Loop over each symmetry operator, applying it to each k-point
        ! Rotate the k-point
        roKpt = matmul(SymOps(:,:,iOp),urKpt)
        write(*,'(/,"Operator:",i3)') iOp 
        write(*,'(3(1x,f7.3))') (SymOps(i,:,iOp),i=1,3)
        write(*,'("roKpt: ",3(f6.3,1x),i3)') roKpt
        ! Make sure the rotated k-point is still part of the kgrid. If not, cycle
        call bring_into_cell(roKpt, InvR, R, eps)
        write(*,'("icKpt: ",3(f6.3,1x),i3)') roKpt
        ! Unshift the k-point temporarily to check if it still is a member of the lattice
        roKpt = roKpt - shift
        if (.not. equal(matmul(invK,roKpt), nint(matmul(invK,roKpt)), eps)) then
           cycle
        endif
        ! Reshift the k-point before finding its index
        roKpt = roKpt + shift
        write(*,'("shKpt: ",3(f6.3,1x),i3)') roKpt
        idx = findKptIndex(roKpt,InvK,L,D)
        write(*,'("Op#",1x,i2,5x,"rkpt: ",3(f6.3,1x),5x,"idx: ",i4)') iOp,roKpt,idx
        ! Is this a new addition to the orbit?
        if (hashTable(idx)==0) then
           hashTable(idx) = cOrbit
           ! if so increase the number of members of this orbit by 1
           iWt(cOrbit)=iWt(cOrbit)+1 
          endif
     enddo
  enddo                                                  
    write(*,*) "Hash table:"
    do i = 1, nUr
       write(*,'("kpt#:",i3,3x,"index:",i3)') i, hashTable(i)
    enddo
    ! Now that we have the hash table populated, make a list of the irreducible k-points
    ! and their corresponding weights.
    sum = 0
    allocate(ReducedList(cOrbit,3),weights(cOrbit))
    weights = iWt(1:cOrbit)
    do i = 1, cOrbit
       sum = sum + weights(i)
       ReducedList(i,:) = UnreducedKpList(iFirst(i),:)
    enddo
    ! Fail safe checks
    if (any(hashTable==0)) then
       write(*,*) "ERROR: (symmetryReduceKpointList in generateKpoints.f90)"
       write(*,*) "At least one k-point in the unreduced list was not included in &
            &one of the orbits of the symmetry group."
       stop
    endif
    if (size(UnreducedKpList,1)/real(size(weights)) > size(SymOps,3)) then
       write(*,*) "ERROR: (symmetryReduceKpointList in generateKpoints.f90)"
       write(*,*) "The ratio of unreduced to reduced kpoints is larger than the size &
            &of the point group."
       stop
    endif
    if (sum /= nUR) then
       write(*,*) "ERROR: (symmetryReduceKpointList in generateKpoints.f90)"
       write(*,*) "The weighted sum of reduced k-points is not equal to the number of&
            & unreduced k-points."
       print*,sum,nuR
       print*,weights
       stop
    endif
  CONTAINS
    ! This function takes a k-point in Cartesian coordinates and "hashes" it into a
    ! single number, corresponding to its place in the k-point list.
    function findKptIndex(kpt, InvK, L, D)
      integer              :: findKptIndex ! index of the k-point (base 10, 1..n)
      ! The k-point, inverse of the k-grid generating vecs
      real(dp), intent(in) :: kpt(3), InvK(3,3) 
      ! Left transform for SNF conversion, diagonal of SNF
      integer,  intent(in) :: L(3,3), D(3) 


      real(dp) :: gpt(3)

      gpt = matmul(InvK,kpt) ! k-point is now in lattice coordinates
!!      if (.not. equal(gpt, nint(gpt), eps)) then ! kpt is not a lattice point of K
!!         write(*,*) "ERROR: (findKptIndex in kpointGeneration)"
!!         write(*,*) "The k-point is not a lattice point of the generating vectors."
!!         stop
!!      END if
      ! Convert the k-point to group coordinates and bring into first unit cell
      gpt = modulo(matmul(L,nint(gpt)),D)
      ! Convert from "group" coordinates (a 3-vector) to single base-10 number
      ! between 1 and nUr
      findKptIndex = gpt(1)*D(2)*D(3) + gpt(2)*D(3) + gpt(3) + 1  ! Hash of the kpt
    END function findKptIndex

  END subroutine symmetryReduceKpointList
  
END MODULE kpointGeneration
