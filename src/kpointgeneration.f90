MODULE kpointgeneration
  use num_types
  use numerical_utilities, only: equal
  use vector_matrix_utilities
  use rational_mathematics, only: HermiteNormalForm, SmithNormalForm_li
  use symmetry
  use group_theory
  implicit none
  private
  public generateIrredKpointList, generateFullKpointList, symmetryReduceKpointList, &
       mapKptsIntoBZ, findQpointsInZone
CONTAINS
  !!<summary> Move a list of k-points into the first Brillioun zone: by applying
  !!translational symmetry, find the set of k-points equivalent to the input set that
  !!is closer to the origin than any other equivalent set. The input set is modified by
  !!this routine. </summary>
  !!<parameter regular="true" name="R"> Matrix of reciprocal lattice vectors as columns of
  !!a 3x3 array. </parameter>
  !!<parameter regular="true" name="KpList"> The list of k-points. </parameter>
  !!<parameter regular="true" name="reps_"> A finite precision
  !!parameter for relative tolerances. (optional) </parameter>
  !!<parameter regular="true" name="aeps_"> A finite precision
  !!parameter for absolute tolerances. (optional) </parameter>
  !!<local name="minkedR"> The basis of the reciprocal lattice in Minkowski space. </local>
  subroutine mapKptsIntoBZ(R, KpList, reps_, aeps_)
    !! <local name="minkedR" dimension="(3,3)"> "The basis of the reciprocal lattice in
    !! Minkowski space. </local>
    real(dp), intent(in)   :: R(3,3)
    real(dp), intent(inout):: KpList(:,:) ! First index over k-points, second coordinates
    real(dp), intent(in), optional :: reps_, aeps_
    real(dp)  :: minkedR(3,3), kpt(3), shift_kpt(3), minkedRinv(3,3), Rinv(3,3), M(3,3)
    real(dp)  :: Minv(3,3)
    real(dp)  :: minLength, length, reps, aeps
    integer   :: ik, nk, i, j, k
    logical   :: err ! flag for catching errors
    
    if(.not. present(reps_)) then
       reps = 1e-10_dp
    else
       reps =  reps_
    endif
    if(.not. present(aeps_)) then
       aeps = 1e-10_dp
    else
       aeps =  aeps_
    endif

    call matrix_inverse(R, Rinv, err, reps)
    if (err) then
       write(*,*) "ERROR: (mapKptsIntoBZ in generateKpoints.f90):"
       write(*,*) "The reciprocal lattice vectors are linearly dependent."
       stop
    endif
    
    call minkowski_reduce_basis(R, minkedR, reps)
    call matrix_inverse(minkedR, minkedRinv,err)
    if (err) then
       write(*,*) "ERROR: (mapKptsIntoBZ in generateKpoints.f90):"
       write(*,*) "The Minkowski reduced lattice vectors are linearly dependent."
       stop
    endif
    
    nk = size(KpList, 1)
    call matrix_inverse(minkedR, minkedRinv, err)
    M = matmul(minkedRinv, R)
    
    if (.not. equal(M, nint(M), reps, aeps)) then
       write(*,*) "ERROR: (mapKptsIntoBZ in generateKpoints.f90):"
       write(*,*) "The Minkowski-reduced basis vectors and the reciprocal& 
            & lattice vectors define different lattices."
       write(*,*) "The Minkowski-reduced basis vectors should be integer& 
            & combinations of the reciprocal lattice vectors."
       stop
    endif
    
    call matrix_inverse(M, Minv, err)
    if (err) then
       write(*,*) "ERROR: (mapKptsIntoBZ in generateKpoints.f90):"
       write(*,*) "The Minkowski transformation matrix is non-invertible."
       stop
    endif
    
    if (.not. equal(Minv, nint(Minv), reps, aeps)) then
       write(*,*) "ERROR: (mapKptsIntoBZ in generateKpoints.f90):"
       write(*,*) "The Minkowski-reduced basis vectors and the reciprocal& 
            & lattice vectors define different lattices."
       write(*,*) "The reciprocal lattice vectors should be integer& 
            & combinations of the Minkowski-reduced basis vectors."
       stop
    endif
    
    do ik = 1, nk
       kpt = KpList(ik,:)
       ! Move the k-point into the first unit cell in the Minkowski basis.
       call bring_into_cell(kpt, minkedRinv, minkedR, aeps)
       KpList(ik,:) = kpt
       minLength = norm(kpt)
       
       ! Look at the eight translationally equivalent k-points in the unit cells that
       ! have a vertex at the origin to see if one of them is closer.
       do i = -1, 0
          do j = -1, 0
             do k = -1, 0
                shift_kpt = i*minkedR(:,1) + j*minkedR(:,2) + k*minkedR(:,3) + kpt
                length = norm(shift_kpt)
                
                if((length + aeps) < minLength) then
                   minLength = length
                   KpList(ik,:) = shift_kpt
                endif
             enddo
          enddo
       enddo
    enddo
  endsubroutine mapKptsIntoBZ
  
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
  !!<parameter regular="true" name="reps_"> A finite precision
  !!parameter for relative tolerances. (optional) </parameter>
  !!<parameter regular="true" name="aeps_"> A finite precision
  !!parameter for absolute tolerances. (optional) </parameter>
  !!<parameter regular="true" name="A">The real space lattice vectors.</parameter>
  !!<parameter regular="true" name="AtBas">The atomic basis in lattice coordinates.</parameter>
  !!<parameter regular="true" name="at">Atomic occupancy list.</parameter>
  !!<parameter regular="true" name="symm_">A flag to control the
  !!symmetries to include. 0 means use all symmetries, 1 means only
  !!use the space group, 2 means only use time reversal symmetry, and
  !!3 means use no symmetry.</parameter>
  !!<parameter name="err_" regular="true">Only present if SNF
  !!overfolws are being checked for, returns 1 in case of
  !!overflow.</parameter>
  subroutine generateIrredKpointList(A,AtBas,at,K, R, kLVshift, IrrKpList, weights, reps_, &
       aeps_, symm_, err_)
    real(dp), intent(in) :: K(3,3), A(3,3)
    real(dp), intent(in) :: R(3,3)
    real(dp), intent(in) :: kLVshift(3)
    real(dp), pointer    :: IrrKpList(:,:)
    real(dp), allocatable :: AtBas(:,:)
    integer, pointer     :: weights(:)
    real(dp), intent(in), optional :: reps_, aeps_
    integer, intent(inout)  :: at(:)
    real(dp), pointer    :: KpList(:,:)
    real(dp), allocatable    :: pgOps(:,:,:), trans(:,:)
    real(dp)             :: reps, aeps
    integer, intent(out), optional :: err_
    integer, intent(in) , optional :: symm_
    integer :: symm

    if(.not. present(reps_)) then
       reps = 1e-10_dp
    else
       reps =  reps_
    endif
    if(.not. present(aeps_)) then
       aeps = 1e-10_dp
    else
       aeps =  aeps_
    endif
    if (.not. present(symm_)) then
       symm = 0
    else
       symm = symm_
    end if

    if (symm == 0) then
       call get_spaceGroup(A, at, AtBas, pgOps, trans, .true., reps)
       call get_fullSpaceGroup(pgOps, reps, aeps)
    else if (symm == 1) then
       call get_spaceGroup(A, at, AtBas, pgOps, trans, .true., reps)
    else if (symm == 2) then
       allocate(pgOps(3,3,2))
       pgOps(:,:,1) = reshape((/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp/), (/3,3/))
       pgOps(:,:,2) = reshape((/-1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -1.0_dp/), (/3,3/))
    else
       allocate(pgOps(3,3,1))
       pgOps(:,:,1) = reshape((/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp/), (/3,3/))
    end if

    call generateFullKpointList(K, R, kLVshift, KpList, reps_=reps, aeps_=aeps)
    if (present(err_)) then
       call symmetryReduceKpointList(K, R, kLVshift, KpList, pgOps, IrrKpList, weights, &
            reps_=reps, aeps_=aeps, err_=err_)
    else 
       call symmetryReduceKpointList(K, R, kLVshift, KpList, pgOps, IrrKpList, weights, &
            reps_=reps, aeps_=aeps)
    end if
    deallocate(KpList,pgOps,trans)
  endsubroutine generateIrredKpointList

  !!<summary>Generates the full space group from the crystal space
  !!group.</summary>
  !!<usage>Pass in 3D array containing the matrices of the spcae
  !!group. The operations will be taken in all possible pairs to
  !!generate new group elements along with the inverse operation. The
  !!array of operations will be added to until the list of
  !!permutations has closure and constitutes a group.</usage>
  !!<parameter name="g" regular="true">On entry contains the
  !!generators, on exit contains the whole group </parameter>
  !!<parameter regular="true" name="reps_"> A finite precision
  !!parameter for relative tolerances. (optional) </parameter>
  !!<parameter regular="true" name="aeps_"> A finite precision
  !!parameter for absolute tolerances. (optional) </parameter>
  subroutine get_fullSpaceGroup(g, reps_, aeps_)
    real(dp), allocatable :: g(:,:,:)
    real(dp), optional :: reps_, aeps_

    real(dp) :: inv(3,3), new_op(3,3)
    real(dp), allocatable :: new_group(:,:,:)
    integer :: i, oG, j, nG
    logical :: inc_inv, new
    real(dp) :: reps, aeps

    if(.not. present(reps_)) then
       reps = 1e-10_dp
    else
       reps =  reps_
    endif
    if(.not. present(aeps_)) then
       aeps = 1e-10_dp
    else
       aeps =  aeps_
    endif

    inv = reshape((/-1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -1.0_dp/), (/3,3/))

    oG = size(g,3)
    ! Check if inversion symmetry is already in the group.
    inc_inv = .False.
    do i = 1, oG
       if (equal(g(:,:,i), inv, reps, aeps)) then
          inc_inv = .True.
          exit
       end if
    end do

    if (.not. inc_inv) then
       !If inversion symmetry isn't in the group then we need to add
       !it and the inverse of each operator in the group, i.e., the
       !group could potentially double in size.
       nG = 0
       allocate(new_group(3,3, oG*2))
       new_group(:,:,1:oG) = g
       do i=1, oG
          new_op = matmul(g(:,:,i), inv)
          new = .True.
          do j=1, oG
             if (equal(new_op, g(:,:,j), reps, aeps)) then
                new = .False.
                exit
             end if
          end do

          ! Verify that the new operator isn't already in the group
          ! just to be safe.
          if (new) then
             nG = nG + 1
             new_group(:,:,oG+nG) = new_op
          end if
       end do

       ! replace the old group with the new one.
       deallocate(g)
       allocate(g(3,3,oG+nG))
       g = new_group
       deallocate(new_group)       
    end if
  end subroutine get_fullSpaceGroup
  
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
  !!<parameter regular="true" name="reps_"> A finite precision
  !!parameter for relative tolerances. (optional) </parameter>
  !!<parameter regular="true" name="aeps_"> A finite precision
  !!parameter for absolute tolerances. (optional) </parameter>
  subroutine generateFullKpointList(K, R, kLVshift, KpList, reps_, aeps_)
    real(dp), intent(in) :: K(3,3)
    real(dp), intent(in) :: R(3,3)
    real(dp), intent(in) :: kLVshift(3)
    real(dp), pointer    :: KpList(:,:)
    real(dp), intent(in), optional :: reps_, aeps_

    real(dp) :: Kinv(3,3) ! Inverse of the k-grid cell 
    real(dp) :: Rinv(3,3) ! Inverse of reciprocal cell
    real(dp) :: reps, aeps ! Finite precision parameter
    real(dp) :: cartShift(3), bicCartShift(3) ! k-grid shifts in Cartesian coordinates
    integer  :: S(3,3), H(3,3), B(3,3) ! Integer matrices for HNF conversion
    integer  :: n ! Number of k-points (i.e., index of kgrid lattice, a.k.a. volume factor)
    integer  :: a, c, f ! Diagonal entries of the HNF matrix
    integer  :: iK, jK, kK ! loop counters of k-point generating vectors
    integer  :: idx ! index (ordinal number) of k-point
    logical  :: err ! flag for catching errors
    real(dp) :: intMat(3,3)
    
    if(.not. present(reps_)) then
       reps = 1e-10_dp
    else
       reps =  reps_
    endif
    if(.not. present(aeps_)) then
       aeps = 1e-10_dp
    else
       aeps =  aeps_
    endif

    ! Check for valid inputs
    if (abs(determinant(K)) > abs(determinant(R))+aeps) then
       write(*,*) "ERROR (generateFullKpointList in generateKpoints.f90):"
       write(*,*) "The k-point generating lattice vectors have a unit cell &
            &larger than the reciprocal lattice."
       write(*,*) abs(determinant(K)), abs(determinant(R))
       stop
    endif
    call matrix_inverse(K,Kinv,err)
    if (err) then
       write(*,*) "ERROR: (generateFullKpointList in generateKpoints.f90):"
       write(*,*) "The kgrid generating vectors are linearly dependent."
       stop
    endif

    intMat = matmul(Kinv,R)
    
    if (.not. equal(matmul(Kinv,R), nint(matmul(Kinv,R)), reps, aeps)) then
       write(*,*) "ERROR: (generateFullKpointList in generateKpoints.f90):"
       stop "The point generating vectors and the reciprocal lattice are incommensurate."
    endif
    
    call matrix_inverse(R,Rinv,err)
    
    if (err) then       
       write(*,*) "ERROR: (generateFullKpointList in generateKpoints.f90):"
       write(*,*) "The reciprocal lattice vectors are linearly dependent."
       stop
    endif

    ! Put the shift in Cartesian coordinates.
    cartShift = matmul(K,kLVshift)

    ! The k-grid shift, but in the unit cell
    bicCartShift = cartShift 
    
    call bring_into_cell(bicCartShift, Kinv, K, aeps)
    
    if (.not. equal(cartShift, bicCartShift, reps, aeps)) then
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
             ! ordinal counter for k-points;
             ! converted to base-10 from mixed-radix number of HNF diagonal entries
             idx = f*c*iK + f*jK + kK + 1
             
             ! compute Cartesian coordinates of the k-point
             KpList(idx,:) = matmul(K,(/iK, jK, kK/)) + cartShift
          enddo
       enddo
    enddo
 
    ! Bring each point in the list into the first BZ.
    ! (Originally, each point was mapped into the first unit cell.)
    call mapKptsIntoBZ(R, KpList, reps_=reps, aeps_=aeps)
    
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
  !!<parameter regular="true" name="reps_"> A finite precision
  !!parameter for relative tolerances. (optional) </parameter>
  !!<parameter regular="true" name="aeps_"> A finite precision
  !!parameter for absolute tolerances. (optional) </parameter>
  !!<parameter name="err_" regular="true">Only present if SNF
  !!overfolws are being checked for, returns 1 in case of
  !!overflow.</parameter>
  subroutine symmetryReduceKpointList(K, R, kLVshift, UnreducedKpList, SymOps, &
       ReducedList, weights, reps_, aeps_, err_)
    real(dp), intent(in) :: K(3,3), R(3,3) 
    real(dp), intent(in) :: kLVshift(3) 
    real(dp), intent(in) :: UnreducedKpList(:,:) 
    real(dp), intent(in) :: SymOps(:,:,:) ! Last slot is op # index, first two are 3x3
    real(dp), pointer    :: ReducedList(:,:)
    integer, pointer     :: weights(:) 
    real(dp), optional   :: reps_, aeps_
    integer, intent(out), optional :: err_

    integer :: iOp, nOps, iUnRdKpt, nUR, cOrbit, idx, i, sum
    integer :: hashTable(size(UnreducedKpList,1))
    integer :: iFirst(size(UnreducedKpList,1)), iWt(size(UnreducedKpList,1))
    real(dp):: InvK(3,3) ! Inverse of kgrid matrix
    real(dp):: InvR(3,3) ! Inverse of reciprocal lattice
    real(dp):: urKpt(3), roKpt(3) ! unrotated k-point, rotated k-point
    real(dp):: shift(3) ! Shift of k-grid lattice in Cartesian coordinates
    integer :: N(3,3) ! Integer transformation that takes K to R

    ! HNF, SNF transform matrices, SNF, diag(SNF)
    integer :: H(3,3), S(3,3), B(3,3)
    integer(li) :: L(3,3), D(3), Ri(3,3)
    real(dp):: reps, aeps
    logical :: err
    integer :: zz
    integer :: intSymOp(3,3,size(SymOps,3))
    
    if(.not. present(reps_)) then
       reps = 1e-10_dp
    else
       reps =  reps_
    endif

    if(.not. present(aeps_)) then
       aeps = 1e-10_dp
    else
       aeps =  aeps_
    endif

    call matrix_inverse(K, InvK, err, aeps)

    if (err) then
       write(*,*) "ERROR: (symmetryReduceKpointList in generateKpoints.f90)"
       stop "The k-grid vectors that were passed in are linearly dependent."
    END if

    call matrix_inverse(R, InvR, err, aeps)
    if (err) then
       write(*,*) "ERROR: (symmetryReduceKpointList in generateKpoints.f90)"
       stop "The reciprocal lattice vectors are linearly dependent."
    END if
    
    ! Make sure kgrid is commensurate with reciprocal cell
    if (.not. equal(matmul(InvK,R), nint(matmul(InvK,R)), reps, aeps)) then
       write(*,*) "ERROR: (symmetryReduceKpointList in generateKpoints.f90):"
       write(*,*) "The point generating vectors and the reciprocal lattice are &
            &incommensurate."
       write(*,*) "The reciprocal lattice should be an integer multiple of&
            & the generating lattice vectors."
       stop
    endif
        
    ! Make sure that kgrid is no bigger than the reciprocal cell.
    ! This determinant check could cause problems for 1-kpoint cases because
    ! eps isn't scaled by the size of the elements the matrices K or R.
    If (abs(determinant(K)) > abs(determinant(R))+aeps) then
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
    call HermiteNormalForm(N,H,B)
    
    ! Left side of transform will be used later, right side (Ri) will not be.
    if (present(err_)) then
       call SmithNormalForm_li(H,L,S,Ri,err_=err_)
       if (err_ == 1) then
          return
       end if
    else
       call SmithNormalForm_li(H,L,S,Ri)       
    end if

    ! Diagonal of the SNF.
    D = (/S(1,1),S(2,2),S(3,3)/)

    ! Count how many unique orbits there are. I.e., the number of unique kpoints
    cOrbit = 0
    
    ! Number of symmetry operators in the point group
    nOps = size(SymOps,3)
    
    ! Number of unreduced kpoints
    nUR = size(UnreducedKpList,1)

    ! Use a hash table to keep track of symmetrically-equivalent kpoints.
    hashTable = 0

    ! Store the index of the first kpoint in each orbit.
    iFirst = 0

    ! Count the number of members of each orbit.
    iWt = 0
    ! Convert the rotation matrix to integer form
    do zz=1,size(SymOps,3)
       intSymOp(:,:,zz) = nint(matmul(InvR, matmul(SymOps(:,:,zz), R)))
    end do

    do iUnRdKpt = 1,nUR ! Loop over each k-point and mark off its symmetry brothers
       zz = 0       
       urKpt = UnreducedKpList(iUnRdKpt,:) ! unrotated k-point (shorter name for clarity)
       idx = findKptIndex(urKpt-shift, InvK, L, D, rtol_=reps, atol_=aeps)
       
       if (hashTable(idx)/=0) cycle ! This k-point is already on an orbit, skip it

       cOrbit = cOrbit + 1
       hashTable(idx) = cOrbit
       iFirst(cOrbit) = iUnRdKpt
       iWt(cOrbit) = 1
       
       do iOp = 1,nOps ! Loop over each symmetry operator, applying it to each k-point
          
          ! Rotate the k-point
          roKpt = matmul(R, matmul(intSymOp(:,:,iOp), matmul(InvR, urKpt)))
          
          ! Make sure the rotated k-point is still part of the kgrid. If not, cycle
          call bring_into_cell(roKpt, InvR, R, aeps)
          
          ! Unshift the k-point to check if it still is a member of the lattice. The
          ! reason for removing the shift is our indexing method doesn't work if the
          ! grid is moved off the origin.
          roKpt = roKpt - shift
          if (.not. equal(matmul(invK,roKpt), nint(matmul(invK,roKpt)), reps, aeps)) then
             cycle
          endif

          ! Find the index of this k-point.
          idx = findKptIndex(roKpt, InvK, L, D, rtol_=reps, atol_=aeps)
          
          ! Is this a new addition to the orbit?
          if (hashTable(idx)==0) then
             
             hashTable(idx) = cOrbit
             
             ! if so increase the number of members of this orbit by 1
             iWt(cOrbit)=iWt(cOrbit)+1
          endif
       enddo
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
          ! stop
       endif
    enddo

    ! Check for closure on the orbits
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
    function findKptIndex(kpt, InvK, L, D, rtol_, atol_)
      
      ! Index of the k-point (base 10, 1..n)
      integer              :: findKptIndex
      
      ! The k-point, matrix inverse of the k-grid generating vecs
      real(dp), intent(in) :: kpt(3), InvK(3,3)
      
      ! Left transform for SNF conversion, diagonal of SNF
      integer(li),  intent(in) :: L(3,3), D(3)
      
      ! Rounding parameter that determines the amount of memory used
      ! to store the number being rounded.
      real(dp), intent(in), optional  :: rtol_, atol_
      
      ! Relative and absolute tolerances
      real(dp) :: rtol, atol
      
      ! The k-point in lattice (grid) coordinates.
      real(dp) :: gpt(3)
      integer(li) ::int_gpt(3)

      ! Set the relative tolerance.
      if(.not. present(rtol_)) then
         rtol = 1e-5_dp
      else
         rtol =  rtol_
      endif
      
      ! Set the absolute tolerance.
      if(.not. present(atol_)) then
         atol = 1e-6_dp
      else
         atol =  atol_
      endif      
      
      ! Put the k-point in lattice (grid) coordinates.
      gpt = matmul(InvK, kpt)
      
      ! Make sure the k-point is a lattice point of K. 
      ! The value of eps, 1E-3, here was empirically chosen when the code was initially tested.
      ! If eps is too tight, this check often fail for very fine meshes. 
      if (.not. equal(gpt, int(nint(gpt,dp)), 1E-3_dp, atol)) then
         write(*,*) "ERROR: (findKptIndex in kpointGeneration)"
         write(*,*) "The k-point is not a lattice point of the generating vectors."
         stop
      END if

      ! Since this k-point is part of the grid, within some tolerance (which was verified
      ! in the previous check), we can safely remove any finite precision errors in the
      ! calculation of gpt and thereby avoid error propagation.
      int_gpt = nint(gpt,dp)
            
      ! Convert the k-point to group coordinates and bring into the first unit cell.      
      gpt = modulo(matmul(L, int_gpt), D)
      
      ! Convert from "group" coordinates (a 3-vector) to single base-10 number
      ! between 1 and nUr      
      findKptIndex = int(gpt(1)*D(2)*D(3) + gpt(2)*D(3) + gpt(3) + 1, 4)  ! Hash of the kpt
      
    END function findKptIndex
    
  END subroutine symmetryReduceKpointList

  !!<summary> Find which Q points in a phonon calculation are
  !! specified ("pinned down") by the supercell
  !! calculation. </summary>
  !!<parameter regular="true" name="Avecs"> Lattice vectors of the
  !!primitive cell </parameter>
  !!<parameter regular="true" name="Bvecs"> Lattice vectors of the
  !!supercell</parameter>
  !!<parameter regular="true" name="n"> The number of q points in the
  !!reciprocal cell of Avecs. This is det(B)/det(A) but we need to
  !!know it ahead of time if we are going to f90-wrap this routine. It
  !!corresponds to the number of qpoints that we will return, so we
  !!need it to avoid an intent(out) variable with an unknown
  !!size. </parameter>
  !!<parameter regular="true" name="qPoints"> (output) The list of Q points </parameter>
  !!<parameter regular="true" name="reps_"> A finite precision
  !!parameter for relative tolerances. (optional) </parameter>
  !!<parameter regular="true" name="aeps_"> A finite precision
  !!parameter for absolute tolerances. (optional) </parameter>
  SUBROUTINE findQpointsInZone(Avecs, Bvecs, n, qPoints, reps_, aeps_)
    real(dp), intent(in), dimension(3,3) :: Avecs, Bvecs
    integer, intent(in)                  :: n
    real(dp), intent(inout)              :: qPoints(n, 3)
    real(dp), intent(in), optional       :: reps_, aeps_


    real(dp) :: reps, aeps, shift(3)
    real(dp), dimension(3,3)             :: Ainv, Binv
    real(dp), pointer :: qList(:,:)
    
    if(.not. present(reps_)) then
       reps = 1e-10_dp
    else
       reps =  reps_
    endif
    if(.not. present(aeps_)) then
       aeps = 1e-10_dp
    else
       aeps =  aeps_
    endif
    
    if (.not. equal(determinant(Bvecs)/determinant(Avecs), n, reps, aeps)) then
       stop "Number of q points expected, 'n', doesn't match the relative sizes of the two lattices."
    end if
    shift= (/0._dp, 0._dp, 0._dp/)
    call matrix_inverse(Avecs, Ainv)
    call matrix_inverse(Bvecs, Binv)
    call generateFullKpointList(Binv, Ainv, shift, qList, reps_=reps, aeps_=aeps)
    
    if (size(qList,1)/= n) stop "Found the wrong number of q points (or 'n' was wrong)"
    qPoints = qList
    call mapKptsIntoBZ(Ainv, qpoints, reps_=reps, aeps_=aeps)
  END SUBROUTINE findQpointsInZone
  
END MODULE kpointgeneration
