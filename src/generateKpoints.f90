  MODULE kpointGeneration
  use num_types
  use numerical_utilities, only: equal
  use vector_matrix_utilities
  use rational_mathematics, only: HermiteNormalForm, SmithNormalForm
  use symmetry
  implicit none

  private
  public generateIrredKpointList, generateFullKpointList, symmetryReduceKpointList, &
       & mapKptsIntoFirstBZ
CONTAINS
  !!<summary> Move a list of k-points into the first Brillioun zone: by applying
  !!translational symmetry, find the set of k-points equivalent to the input set that
  !!is closer to the origin than any other equivalent set. The input set is modified by
  !!this routine. </summary>
  !!<parameter regular="true" name="R"> Matrix of reciprocal lattice vectors as columns of
  !!a 3x3 array. </parameter>
  !!<parameter regular="true" name="KpList"> The list of k-points. </parameter>
  !!<parameter regular="true" name="eps_"> A finite precision parameter. (optional)
  !!</parameter>
  !!<local name="minkedR"> The basis of the reciprocal lattice in Minkowski space. </local>
  subroutine mapKptsIntoFirstBZ(R, KpList, eps_)
    !! <local name="minkedR" dimension="(3,3)"> "The basis of the reciprocal lattice in Minkowski space. </local>
    Real(dp), intent(in)   :: R(3,3)
    real(dp), intent(inout):: KpList(:,:) ! First index over k-points, second coordinates
    real(dp), intent(in), optional :: eps_
    real(dp)  :: minkedR(3,3), kpt(3), shift_kpt(3), minkedRinv(3,3), Rinv(3,3), M(3,3)
    real(dp)  :: Minv(3,3)
    real(dp)  :: minLength, length, eps
    integer   :: ik, nk, i, j, k
    logical   :: err ! flag for catching errors
    
    if(.not. present(eps_)) then
       eps = 1e-10_dp
    else
       eps =  eps_
    endif
 
    call matrix_inverse(R, Rinv, err)
    if (err) then
       write(*,*) "ERROR: (mapKptsIntoFirstBZ in generateKpoints.f90):"
       write(*,*) "The reciprocal lattice vectors are linearly dependent."
       stop
    endif

    call minkowski_reduce_basis(R, minkedR, eps)
    call matrix_inverse(minkedR, minkedRinv, err)
    if (err) then
       write(*,*) "ERROR: (mapKptsIntoFirstBZ in generateKpoints.f90):"
       write(*,*) "The Minkowski reduced lattice vectors are linearly dependent."
       stop
    endif
    
    nk = size(KpList, 1)
    call matrix_inverse(minkedR, minkedRinv, err)
    M = matmul(minkedRinv, R)
    if (.not. equal(M, nint(M), eps)) then
       write(*,*) "ERROR: (mapKptsIntoFirstBZ in generateKpoints.f90):"
       write(*,*) "The Minkowski-reduced basis vectors and the reciprocal& 
            & lattice vectors define different lattices."
       write(*,*) "The Minkowski-reduced basis vectors should be integer& 
            & combinations of the reciprocal lattice vectors."
       stop
    endif
    
    call matrix_inverse(M, Minv, err)
    if (err) then
       write(*,*) "ERROR: (mapKptsIntoFirstBZ in generateKpoints.f90):"
       write(*,*) "The Minkowski transformation matrix is non-invertible."
       stop
    endif
    
    if (.not. equal(Minv, nint(Minv), eps)) then
       write(*,*) "ERROR: (mapKptsIntoFirstBZ in generateKpoints.f90):"
       write(*,*) "The Minkowski-reduced basis vectors and the reciprocal& 
            & lattice vectors define different lattices."
       write(*,*) "The reciprocal lattice vectors should be integer& 
            & combinations of the Minkowski-reduced basis vectors."
       stop
    endif
    
    do ik = 1, nk
       kpt = KpList(ik,:)
       ! Move the k-point into the first unit cell in the Minkowski basis.
       call bring_into_cell(kpt, minkedRinv, minkedR, eps)
       KpList(ik,:) = kpt
       minLength = norm(kpt)
       
       ! Look at the eight translationally equivalent k-points in the unit cells that
       ! have a vertex at the origin to see if one of them is closer.
       do i = -1, 0
          do j = -1, 0
             do k = -1, 0
                shift_kpt = i*minkedR(:,1) + j*minkedR(:,2) + k*minkedR(:,3) + kpt
                length = norm(shift_kpt)
                ! write(*,'("shift_kpt: ",3(f6.3,1x))') shift_kpt
                if((length + eps) < minLength) then
                   ! write(*,'("n: ",3x,f4.1)') length
                   minLength = length
                   KpList(ik,:) = shift_kpt
                endif
             enddo
          enddo
       enddo
    enddo
  endsubroutine mapKptsIntoFirstBZ
  
  !!<summary>Reduce k-points to irreducible set. Takes a list of translationally distinct,
  !! but not rotationally distinct, k-points and reduces them by the given point group.
  !! </summary>
  !!<parameter regular="true" name="K"> Matrix of grid generating vectors as columns of a
  !! 3x3 array. </parameter>
  !!<parameter regular="true" name="R"> Matrix of reciprocal lattice vectors as columns of
  !! a 3x3 array. </parameter>
  !!<parameter regular="true" name="kLVshift"> Offset (shift) of the k-grid in fractions
  !! of k-grid vectors. </parameter>
  !!<parameter regular="true" name="IrrKpList"> List of symmetry-reduced k-points in
  !! Cartesian coordinates. </parameter>
  !!<parameter regular="true" name="weights"> "Weights" of k-points (length of each orbit).
  !!</parameter>
  !!<parameter regular="true" name="eps_">Finite precision parameter (optional)</parameter>
  !!<parameter regular="true" name="A">The real space lattice vector.</parameter>
  !!<parameter regular="true" name="B_vecs">The basis vectors of the lattice in lattice
  !!coordinates.</parameter>
  !!<parameter regular="true" name="at">Atomic occupancy list.</parameter>
  subroutine generateIrredKpointList(A,B_vecs,at,K, R, kLVshift, IrrKpList, weights, eps_)
    real(dp), intent(in) :: K(3,3), A(3,3)
    real(dp), intent(in) :: R(3,3)
    real(dp), intent(in) :: kLVshift(3)
    real(dp), pointer    :: IrrKpList(:,:), B_vecs(:,:)
    integer, pointer     :: weights(:)
    real(dp), intent(in), optional:: eps_
    integer, intent(inout)  :: at(:)

    real(dp), pointer    :: KpList(:,:)
    real(dp), pointer    :: pgOps(:,:,:), trans(:,:)
    real(dp)             :: eps

    if(.not. present(eps_)) then
       eps = 1e-10_dp
    else
       eps =  eps_
    endif
    
    ! call get_lattice_pointGroup(R, pgOps, eps)
    call get_spaceGroup(A,at,B_vecs,pgOps,trans, .True., eps_=eps)
    call generateFullKpointList(K, R, kLVshift, KpList, eps)
    call symmetryReduceKpointList(K, R, kLVshift, KpList, pgOps, IrrKpList, weights, eps)
  endsubroutine generateIrredKpointList

  !!<summary> Takes two lattices, a generating lattice (K) and the reciprocal lattice (R),
  !! and generates all the points of K that are inside one unit cell of R. </summary>
  !!<parameter name="K" regular="true"> 3x3 matrix. Columns are the generating vectors of
  !! the k-grid </parameter>
  !!<parameter name="R" regular="true"> 3x3 matrix. Columns are the reciprocal lattice
  !! vectors </parameter>
  !!<parameter name="kLVshift" regular="true"> Fractional shift of the k-grid. Given as
  !! fractions of the generating k-grid vectors. </parameter> 
  !!<parameter name="KpList" regular="true"> List of unreduced k-points in Cartesian
  !! coordinates. </parameter>
  subroutine generateFullKpointList(K, R, kLVshift, KpList, eps_)
    real(dp), intent(in) :: K(3,3)
    real(dp), intent(in) :: R(3,3)
    real(dp), intent(in) :: kLVshift(3)
    real(dp), pointer    :: KpList(:,:)
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

    real(dp) :: intMat(3,3)
    integer  :: i

    ! real(dp) :: test1(3,3) ! Inverse of the k-grid cell
    ! real(dp) :: test2(3,3) ! Inverse of the k-grid cell
    ! integer  :: i
    
    if(.not. present(eps_)) then
       eps = 1e-10_dp
    else
       eps =  eps_
    endif

    ! Check for valid inputs
    if (abs(determinant(K)) > abs(determinant(R))+eps) then
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

    ! write(*,'(3("M1: ",3(1x,f50.30),/))') (tmpM1(i,:),i=1,3)
    ! write(*,'(3("M2: ",3(1x,f50.30),/))') (tmpM2(i,:),i=1,3)
    ! write(*, *) equal(matmul(Kinv,R), nint(matmul(Kinv,R)), eps)
    ! write(*,'(3("K: ",3(1x,f11.7),/))') (K(i,:),i=1,3)
    ! test1 = matmul(Kinv,R)
    ! test2 = nint(matmul(Kinv,R))
    ! write(*,'(3("test1: ",3(1x,f11.7),/))') (test1(i,:),i=1,3)
    ! write(*,'(3("test2: ",3(1x,f11.7),/))') (test2(i,:),i=1,3)

    intMat = matmul(Kinv,R)
    write(*,'(3("test1: ",3(1x,f11.7),/))') (intMat(i,:),i=1,3)

    
    if (.not. equal(matmul(Kinv,R), nint(matmul(Kinv,R)), eps)) then
    ! if (any(matmul(Kinv,R) -  nint(matmul(Kinv,R)) > eps)) then
       write(*,*) "ERROR: (generateFullKpointList in generateKpoints.f90):"
       stop "The point generating vectors and the reciprocal lattice are incommensurate."
    endif
    call matrix_inverse(R,Rinv,err)
        if (err) then
       write(*,*) "ERROR: (generateFullKpointList in generateKpoints.f90):"
       write(*,*) "The reciprocal lattice vectors are linearly dependent."
       stop
    endif
    cartShift = matmul(K,kLVshift)
    bicCartShift = cartShift ! The k-grid shift, but bring into cell
    call bring_into_cell(bicCartShift, Kinv, K, eps)
    
    if (.not. equal(cartShift, bicCartShift, eps)) then
       write(*,*) "WARNING: (generateFullKpointList in generateKpoints.f90)"
       write(*,*) "The given shift was outside the first k-grid cell. By translational"
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
             KpList(idx,:) = matmul(K,(/iK, jK, kK/)) + cartShift
          enddo
       enddo
    enddo
    
    ! Bring each k-point into the first unit cell
    do iKP = 1,n
       call bring_into_cell(KpList(iKP,:),Rinv,R,eps)
    enddo

  END subroutine generateFullKpointList

  !!<summary> Reduce k-points to irreducible set. Takes a list of translationally distinct,
  !! but not rotationally distinct, k-points and reduces them by the given point group.
  !!</summary>
  !!<parameter regular="true" name="K"> Matrix of grid generating vectors as columns of a
  !! 3x3 array. </parameter>
  !!<parameter regular="true" name="R"> Matrix of reciprocal lattice vectors as columns of
  !! a 3x3 array. </parameter>
  !!<parameter regular="true" name="kLVshift"> Offset (shift) of the k-grid in fractions
  !! of k-grid vectors. </parameter>
  !!<parameter regular="true" name="UnreducedKpList"> Unreduced k-point list in Cartesian
  !! coordinates. </parameter>
  !!<parameter regular="true" name="SymOps"> Integer rotation matrices in the coordinates
  !! of the lattice basis). </parameter>
  !!<parameter regular="true" name="ReducedList"> List of symmetry-reduced k-points in
  !! Cartesian coordinates. </parameter>
  !!<parameter regular="true" name="weights"> "Weights" of k-points (length of each orbit).
  !! </parameter>
  !!<parameter regular="true" name="eps_"> Finite precision parameter (optional).
  !! </parameter>
  subroutine symmetryReduceKpointList(K, R, kLVshift, UnreducedKpList, SymOps, &
       ReducedList, weights, eps_)
    real(dp), intent(in) :: K(3,3), R(3,3) 
    real(dp), intent(in) :: kLVshift(3) 
    real(dp), intent(in) :: UnreducedKpList(:,:) 
    real(dp), intent(in) :: SymOps(:,:,:) ! Last slot is op # index, first two are 3x3
    real(dp), pointer    :: ReducedList(:,:)
    integer, pointer     :: weights(:) 
    real(dp), optional   :: eps_

    integer :: iOp, nOps, iUnRdKpt, nUR, cOrbit, idx, i, sum
    integer :: hashTable(size(UnreducedKpList,1))
    integer :: iFirst(size(UnreducedKpList,1)), iWt(size(UnreducedKpList,1))
    real(dp):: InvK(3,3) ! Inverse of kgrid matrix
    real(dp):: InvR(3,3) ! Inverse of reciprocal lattice
    real(dp):: urKpt(3), roKpt(3) ! unrotated k-point, rotated k-point
    real(dp):: shift(3) ! Shift of k-grid lattice in Cartesian coordinates
    integer :: N(3,3) ! Integer transformation that takes K to R
    ! HNF, SNF transform matrices, SNF, diag(SNF)
    integer :: H(3,3), L(3,3), Ri(3,3), S(3,3), D(3) 
    real(dp):: eps
    logical :: err
    integer zz
    
    if(.not. present(eps_)) then
       eps = 1e-10_dp
    else
       eps =  eps_
    endif
        
    call matrix_inverse(K, InvK, err, eps)
    if (err) then
       write(*,*) "ERROR: (symmetryReduceKpointList in generateKpoints.f90)"
       stop "The k-grid vectors that were passed in are linearly dependent."
    END if

    call matrix_inverse(R, InvR, err, eps)
    if (err) then
       write(*,*) "ERROR: (symmetryReduceKpointList in generateKpoints.f90)"
       stop "The reciprocal lattice vectors are linearly dependent."
    END if    
    !! Check for valid inputs
    ! Make sure kgrid is commensurate with reciprocal cell
    ! if (any(matmul(InvK,R) -  nint(matmul(InvK,R)) > eps)) then
    if (.not. equal(matmul(InvK,R), nint(matmul(InvK,R)), eps)) then
       write(*,*) "ERROR: (symmetryReduceKpointList in generateKpoints.f90):"
       write(*,*) "The point generating vectors and the reciprocal lattice are &
            &incommensurate."
       write(*,*) "The reciprocal lattice should be an integer multiple of&
            & the generating lattice vectors."
       stop
    endif
        
    ! Make sure that kgrid is no bigger than reciprocal cell
    ! This determinant check could cause problems for 1-kpoint cases because
    ! eps isn't scaled by the size of the elements the matrices K or R.
    If (abs(determinant(K)) > abs(determinant(R))+eps) then
       write(*,*) "ERROR (symmetryReduceKpointList in generateKpoints.f90):"
       write(*,*) "The k-point generating lattice vectors define a unit cell larger&
            & than the reciprocal lattice."
       write(*,*) "This doesn't make sense. There should be at least 1 k-point per &
            &reciprocal cell."
       stop
    endif
    
    ! Put the shift in Cartesian coordinates.
    shift = matmul(K,kLVshift)
    ! Integer transformation matrix that takes K to R, not necessarily HNF at the outset
    N = nint(matmul(InvK,R))
    ! Find the HNF of N, store it in H (B is the transformation matrix)
    call HermiteNormalForm(N,H,L)
    ! Left side of transform will be used later, right side (Ri) will not be.
    call SmithNormalForm(H,L,S,Ri)

    ! Diagonal of the SNF.
    D = (/S(1,1),S(2,2),S(3,3)/)

    ! write(*,'(3(1x,i3))') (H(i,:),i=1,3)
    ! write(*,*) 
    ! write(*,'(3(1x,i3))') (S(i,:),i=1,3)
    ! write(*,'("InvK")')
    ! write(*,'(3(1x,f7.3))') (InvK(i,:),i=1,3)
    cOrbit = 0 ! Count how many unique orbits there are. I.e., the number of unique kpoints
    nOps = size(SymOps,3) ! Number of symmetry operators in the point group
    nUR = size(UnreducedKpList,1) ! Number of unreduced kpoints
    hashTable = 0 ! Use a hash table to keep track of symmetrically-equivalent kpoints
    iFirst = 0 ! store the index of the first kpoint in an orbit
    iWt = 0 ! count the number of members of each orbit
    do iUnRdKpt = 1,nUR ! Loop over each k-point and mark off its symmetry brothers
       zz = 0       
       urKpt = UnreducedKpList(iUnRdKpt,:) ! unrotated k-point (shorter name for clarity)
       idx = findKptIndex(urKpt-shift, InvK, L, D)       
       if (hashTable(idx)/=0) cycle ! This k-point is already on an orbit, skip it
       cOrbit = cOrbit + 1
       hashTable(idx) = cOrbit
       iFirst(cOrbit) = iUnRdKpt
       iWt(cOrbit) = 1
       ! write(*,'("urKpt: ",3(f6.3,1x),i3)') urKpt
       ! write(*,'("index: ",i7)') idx
       do iOp = 1,nOps ! Loop over each symmetry operator, applying it to each k-point
          ! Rotate the k-point
          roKpt = matmul(SymOps(:,:,iOp),urKpt)
          ! write(*,'(/,"Operator:",i3)') iOp
          ! write(*,'(3(1x,f7.3))') (SymOps(i,:,iOp),i=1,3)
          ! write(*,'("shift: ",3(f6.3,1x),i3)') shift
          ! write(*,'("kpoint in lattice coords: ",3(f6.3,1x))') matmul(InvR,roKpt)
          
          ! Make sure the rotated k-point is still part of the kgrid. If not, cycle
          call bring_into_cell(roKpt, InvR, R, eps)
          
          ! write(*,'("roKpt: ",3(f6.3,1x),i3)') roKpt
          ! write(*,'("into cell: ",3(f6.3,1x),i3)') roKpt
          ! write(*,'("icKpt: ",3(f6.3,1x),i3)') roKpt
          ! Unshift the k-point temporarily to check if it still is a member of the lattice
          roKpt = roKpt - shift
          ! write(*,'("no shift: ",3(f6.3,1x),i3)') roKpt          
          
          if (.not. equal(matmul(invK,roKpt), nint(matmul(invK,roKpt)), eps)) then
             zz = zz + 1
             cycle
          endif
          ! write(*,'(3(1x,f7.3))') (SymOps(i,:,iOp),i=1,3)          
          
          idx = findKptIndex(roKpt, InvK, L, D)
          ! write(*,'("rotated index: ",i7)') idx          
          ! Reshift the k-point before finding its index
          roKpt = roKpt + shift
          ! write(*,'("with shift: ",3(f6.3,1x),i3)') roKpt
          
          ! write(*,'("shKpt: ",3(f6.3,1x),i3)') roKpt
          ! idx = findKptIndex(roKpt, InvK, L, D)
          ! write(*,'("Op#",1x,i2,5x,"rkpt: ",3(f6.3,1x),5x,"idx: ",i4)') iOp,roKpt,idx
          ! Is this a new addition to the orbit?
          if (hashTable(idx)==0) then
             ! write(*,'(/,"Operator:",i3)') iOp 
             ! write(*,'(3(1x,f7.3))') (SymOps(i,:,iOp),i=1,3)
             ! write(*,'("roKpt: ",3(f6.3,1x),i3)') roKpt
             ! write(*,'("Op#",1x,i2,5x,"rkpt: ",3(f6.3,1x),5x,"idx: ",i4)') iOp,roKpt,idx
             
             ! print *,"point added"
             hashTable(idx) = cOrbit
             ! if so increase the number of members of this orbit by 1
             iWt(cOrbit)=iWt(cOrbit)+1
          endif
       enddo
       ! print *,'zz ',zz
       ! write(*,'(/,"iWt: ",i4)') iWt(cOrbit)
    enddo
    ! write(*,*) "Hash table:"
    do i = 1, nUr
       ! write(*,'("kpt#:",i3,3x,"index:",i3)') i, hashTable(i)
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

    
    ! ******** Consistency checks ********   
    ! Check that the orbit length divides the group order
    do i = 1, cOrbit
       if (mod(nOps, weights(i)) /= 0) then
          write(*,*) "ERROR: (symmetryReduceKpointList in generateKpoints.f90)"
          write(*,*) "The length of an orbit did not divide the group size"
          write(*,'("Group size:",i3,"   Orbit length:",i3)') nOps, weights(i)
          stop
       endif
    enddo

    ! Check for closure on the orbits
    !
    ! Code still pending...
    
    ! Fail safe checks
    if (any(hashTable==0)) then
       write(*,*) "ERROR: (symmetryReduceKpointList in generateKpoints.f90)"
       write(*,*) "At least one k-point in the unreduced list was not included in &
            &any of the orbits of the symmetry group."
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
       write(*,'("Sum of irr kpoints: ",i7," # of Unred. points:",i7)') sum, nUR
       stop
    endif
    
  CONTAINS
    ! This function takes a k-point in Cartesian coordinates and "hashes" it into a
    ! single number, corresponding to its place in the k-point list.
    function findKptIndex(kpt, InvK, L, D)
      integer              :: findKptIndex ! index of the k-point (base 10, 1..n)
      ! The k-point, matrix inverse of the k-grid generating vecs
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
      ! write(*,*) "index inside", int(gpt(1)*D(2)*D(3) + gpt(2)*D(3) + gpt(3) + 1)
      
      findKptIndex = int(gpt(1)*D(2)*D(3) + gpt(2)*D(3) + gpt(3) + 1)  ! Hash of the kpt
    END function findKptIndex

  END subroutine symmetryReduceKpointList
  
END MODULE kpointGeneration
