0.6.0 (GLWH)
* Added the functionality to pass in a superlattice (not
  crystal) and generate the unreduced kpoint list. The utility is that
  this allows us to see which "q-points" in phonon calculation will be
  calculated explicitly due to the supercell size and shape.
* The function added is currently at the bottom of the file and is
  called "FindQpointsInZone". The function returns the entire list of
  points but moves them into the first BZ.

0.5.12 (JJ)
* Removed print statements from code that were used for testing.
* Since all we care about is the index, removed line that added 
  shift back to the k-point in `symmetryReduceKpointList`.


0.5.11 (GLWH)
* At the end of the generateFullKpointList routine, I
removed the loop that knocked all the kpoints into the first unit
cell. This seem pointless. Instead, there are now mapped into the
first BZ

0.5.10 (JJ)
* The previous fix wasn't really a fix. In the calculation of the
  k-point's unique index associated with its position on the grid,
  errors were being introduced from `L` and `D` not having enough
  digits. These were converted to long integers. Also the k-point in
  grid coordinates had decimal portions removed to avoid error
  propagation.

0.5.9 (JJ)
* Fixed a rounding bug when converting a k-point from Cartesian to
  lattice coordinates.

0.5.8 (WSM)
* Changed the point group finding to finding the space group.

0.5.7 (JJ)
* Updated organization of unittests so that fortpy would be able to
  compile them.

0.5.6 (GLWH)
* Added a check to make sure that the length of each "orbit" of
kpoints divided the order of the point group. 

0.5.5 (JJ)
* Made doc strings clearer by specifying the coordinate system of the
  k-points and clarifying that the lattice vector or grid generating
  vectors are the columns of matrices.
* Added a unit test for 'mapKptsIntoFirstBZ' for body-centered_cubic.

0.5.4 (WSM)
* Removed fortpy.f90 from the repository and from the makefile (since
  it's only needed for unit testing it doesn't make sense to keep it
  in the repo permanently).

0.5.3 (KL)
* Added a unit test for 'mapKptsIntoFirstBZ' for body-centered_cubic.

0.5.2 (JJ)
* Fixed `mapKptsIntoFirstBZ`. The function now takes a list of k-points, maps them to the
  first unit cell in the Minkowski basis, then looks at the translationally equivalent
  k-points in the eight unit cells that have a vertex at the origin to see if any them lie
  closer to the origin.
* Removed previous unit tests and added a few for body-centered cubic.

0.5.1 (K.L.)
* Changed the driver so that it correctly maps the shift into the Minkowski unit cell.
* Added the first simple cubic unit test for the routine 'mapKptsIntoFirstBZ'.

0.5.0 (K.L.)
* Added a new routine that maps a point into the minkowski unit cell. Mapping the grid into
  the Minkowski unit cell is required before symmetry reduction and moving k-points into
  the first Brillouin zone since different points in the grid are in the same orbitals
  after the grid is mapped into the Minkowski unit cell than were in the same orbitals in
  the original grid in the reciprocal unit cell.
* Rearranged the order of steps in the driver. The driver now maps the k-points into the
  Minkowski unit cell, finds the point group of the Minkowski-reduced basis, symmetry
  reduces the grid, and then moves the k-ponts into the first Brillouin zone. This yields
  an accurate list of k-point orbitals in the Brillouin zone.
* Changed the routine 'mapKptsIntoFirstBZ' so that it takes the k-points into the unit cell
  before mapping the grid into Minkowski space. This fixed some issues with incorrect
  mapping.
* Removed the simple cubic and body-centered cubic unit tests for the routine
  'mapKptsIntoFirstBZ' since fixing the routine made these unit tests void.
* Added error checks in 'mapKptsIntoFirstBZ' to ensure that the Minkowski-reduced basis
  vectors and the reciprocal lattice vectors define equivalent lattices.
* Added error check in 'mapKptsIntoFirstBZ' to ensure that the k-points are mapped into the
  Minkowski unit cell.

0.4.3 (JJ)
* Moved the code inside `mapKptsIntoFirstBZ` that determines how many lattice points
  to look in each direction outside the main for loop since it has no k-point
  dependence.
* Made many changes to `mapKptsIntoFirstBZ`: the k-points are placed in the unit cell,
  transformed to Minkowski space, and then moved to the first Brillouin zone in
  Minkowski space.
* The portion of the code in `mapKptsIntoFirstBZ` that searches for the translationally
  equivalent point nearest the origin was removed. The search is always performed within
  the 8 unit cells that surround the origin.
* All that checks in generateKpoints.f90 that verified two lattices were commensurate
  were made consistent by using the function `equal`.
* Added unit tests for the `mapKptsIntoFirstBZ` for a simple cubic lattice.

0.4.2 (K.L.)
* The body-centered_cubic unit tests folder had somehow been removed, so I readded it from
  a previous commit.

0.4.1 (JJ)
* Mapped points into the unit cell before mapping into the first Brillouin zone in
  `mapKptsIntoFirstBZ`.
* The lattice vectors and not the Minkowski reduced lattice vectors were being implemented
  in `mapKptsIntoFirstBZ`. Replaced `R` with `minkedR` in this routine.
* Fixed descriptions of unit tests in generateKpoints.xml.
* It seemed redundant, for example, to have a unit test file named simple_cubic_kpts.in1
  inside the folder simple_cubic. I renamed the files
  
  simple_cubic/unreduced_klist.in.* -> simple_cubic/klist.in.*
  
  simple_cubic/unreduced_klist.in.vasp.* -> simple_cubic/klist.in.vasp*
  
  simple_cubic/simple_cubic_kpts.out.* -> simple_cubic/klist.out.*
  
  simple_cubic/simple_cubic_kpts.out.vasp.* -> simple_cubic/klist.out.vasp*

  simple_cubic/simple_cubic_wts.out.* -> simple_cubic/weights.out.*  
  
  simple_cubic/simple_cubic_wts.out.vasp* -> simple_cubic/weights.in.vasp*
  
  This was applied to all the Bravais lattices.
  
0.4.0 (GLWH)
* Added a new routine that maps a list of k-points (reduced or not) into the first BZ.
  Needs to be unit tested.
  (Jeremy suggests mapping any k-point into the first unit cell before doing check.)

0.3.12 (K.L.)
* Added unit tests that compare k-point reduction to the results obtained from VASP for the
  rhombohedral and triclinic lattices.

0.3.11 (K.L.)
* Made epsilon bigger in the driver because epsilon was not sufficiently large to correctly 
  identify the symmetry operators for hexagonal lattices.
* Added unit tests that compare k-point reduction to that obtained in VASP for the hexagonal
  crystal class.

0.3.10 (JJ)
* Added unit tests that compare k-point reduction to that obtain in VASP for the tetragonal
  crystal classes.
  
0.3.9 (JJ)
* Had to redo the orthorhombic unit tests since the lattice vectors in VASP are on the rows
  instead of the columns. Added unit tests that compare k-point reduction to that obtain in
  VASP for the monoclinic crystal classes.

0.3.8 (JJ)
* Added unit tests that compare k-point reduction to that obtain in VASP for primitive,
  base-centered, body-centered, and face-centered, orthorhombic lattices.

0.3.7 (JJ)
* Added unit tests that compare k-point reduction to that obtain in VASP for simple,
  body-centered, and face-centered, cubic lattices.

0.3.6 (K.L.)
* Added more unit tests for tetragonal lattices.

0.3.5 (K.L.)
* Removed nirrkpts from 'driver.f90'.
* Added unit tests for tetragonal lattices.
* The error for the grid generating vectors being larger than the reciprocal unit cell
  occurred when it shouldn't have when the determinant of the matrix K was negative. This
  same error also didn't occur when it should have when the determinant of the matrix R was
  negative. I fixed this by taking the absolute value of each determinant when this error
  was tested for in each subroutine of 'generatekpoints.f90'.

0.3.4
* Reverted subroutine 'symmetryReduceKpointList' to Revision 0.3.2 (removed optional
  output 'nirrkpts'). Unit tests that compare symmetry reduction to VASP will manually
  be verified to ensure the correct reduction before creating the irreducible k-points to
  compare.

0.3.3
* Added an optional output 'nirrkpts' for 'symmetryReduceKpointList' that gives the
  number of irreducible k-points. The motivation was to be able to create unit tests
  that would compare our k-point reduction to that obtained in VASP.
* Added unit tests that compare k-point reduction to that obtain in VASP for simple
  cubic and face-centered cubic lattices.

0.3.2
* Added 10 unit tests for monoclinic, base-centered monoclinic, hexagonal, rhombohedral,
  and triclinic lattices.

0.3.1
* Added an additional routine, one that combines the "generateFullKpointList" and
"symmetryReduceKpointList" into a single call. This is the use case most users would want.

0.2.9
* Added 10 unit tests for base-centered orthorhombic, body-centered orthorhombic, and
  face-centered orthorhombic.

0.2.8
* Added unit tests for body-centered tetragonal and orthorhombic lattices.

0.2.7
* Added unit tests for face-centered cubic lattices.
* Fixed small bug associated with space between xml documentation and subroutine
  declaration introduced in 0.2.6.

0.2.6
* The warning associated with a shift outside the first k-point cell had the arguments of
  bring_into_cell in the wrong order.
* Changed a few error messages from 'the k-grid vectors are not linearly dependent' to
  'the k-grid vectors are linearly dependent'.
* The k-point index of the unreduced k-point and the rotated k-point were different when
  the rotation operator was the identity. The shift was removed from the unreduced k-point
  whereas it wasn't for the rotated k-point. This was fixed by calculating the rotated
  k-point's index before adding the shift back on.
* Removed repeated code found in symmetryReduceKpointList.
* Fixed an issue where the k-point that represented an orbit was wrong. Instead of having
  iFirst point from cOrbit to the index of the k-point that represented the orbit, changed
  it so that cOrbit pointed to the index of the point in UnreducedKpList.
  ```
  ! iFirst(cOrbit) = idx
  iFirst(cOrbit) = iUnRdKpt
  ```
  
0.2.5
* Enforced line breaks at 90 characters to help readability
* Fixed the O(N^2) problem that Martijn pointed out (fortran intrinsics
  "count" and "minloc" are implicitly O(N^2) ). Actually makes a HUGE difference
  for large cases. The O(N^2) problem was introduced sometime after the initial working
  copy (0.0.2 or so).
* Removed the check for off-kgrid-lattice k-points inside of "findKptIndex"
* Some of the unit tests are not working yet.

0.2.4
* Deleted a portion of the code in symmetryReduceKpointList that was a copy of code in
  findKptIndex. This fixed a bug.
* Started enforcing new lines when line length exceeded 70 characters.
* Tidied up HISTORY.md and findKptIndex.
* Made it so the offset in generateFullKpointList always lies within the first unit cell.
 
* Added a few more unit tests for simple cubic lattices.

0.2.3 
* Added new unit tests for body-centered, cubic lattices.

0.2.2
* Fixed finite precision issues related to eps.
* Added a few unit tests for simple cubic lattices.
* Added an additional fail safe test.

0.2.1
* Bug fix: fixed finite precision error when checking if a rotated kpoint is off the
  k-grid. Was using an epsilon check without using abs() and when the values were
  negative, the check incorrectly failed. Replaced the condition with the "equal" function
  from numerical_utilities module. Should have been using that all along anyway---I coded
  it up so that I wouldn't make these kinds of mistakes! Replaced the line:
 ```if (any((matmul(invK,roKpt) - nint(matmul(invK,roKpt))) > eps)) then```
 with
 ```if (.not. equal(matmul(invK,roKpt), nint(matmul(invK,roKpt)), eps)) then```
      
0.2.0
* Bug fixes, redefinition of shift. Arbitrary shifts seem to do the right thing. If a
  rotation maps a kpoint onto a point not in original list, its skipped (as we want---no
  expansion of the original list). Earlier cases seem to still work with these changes.
  Time for the unit tests now. 

0.1.1
* Included a write-up from Rod that describes how the SNF is used to map points into
  "group coordinates". This is a key idea for the O(N) algorithm, which uses the mapping
  as a hash function for the kpoint list.

0.1.0
* Added functionality for an arbitrary shift. Seems to be working (tried several cases by
  hand). Next thing we need is a large collection of unit tests.

<<< Note, still need to add functionality for skipping rotations that move a k-grid point
off the k-grid lattice. This can happen for strange shifts---which we should avoid, but it
would be nice for the code to do the right thing anyhow. >>> 

0.0.3


* Somehow, and old, incomplete copy got pushed and v. 0.0.2, if it was complete, got lost.
  This version seems to be working again, with the same functionality.


0.0.2
* Added files from laptop. First working copy.

0.0.1

* Initial import of rough, working copy of the folding code. Relies on routines in symlib
  (getPointGroup, SNF, HNF, etc.)

*Somehow, an old, incomplete copy got pushed and v. 0.0.2, if it was complete, got lost.
 This version seems to be working again, with the same functionality.

