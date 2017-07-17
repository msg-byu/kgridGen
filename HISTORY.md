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
<<<<<<< HEAD
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
=======
* Bug fix: fixed finite precision error when checking if a rotated
  kpoint is off the k-grid. Was using an epsilon check without
  using abs() and when the values were negative, the check
  incorrectly failed. Replaced the condition with the "equal"
  function from numerical_utilities module. Should have been
  using that all along anyway---I coded it up so that I wouldn't
  make these kinds of mistakes!
  Replaced the line:
  ```
  if (any((matmul(invK,roKpt) - nint(matmul(invK,roKpt))) > eps)) then
  ```
  with
  ```
  if (.not. equal(matmul(invK,roKpt), nint(matmul(invK,roKpt)), eps)) then
  ```
      
0.2.0
* Bug fixes, redefinition of shift. Arbitrary shifts seem to do
  the right thing. If a rotation maps a kpoint onto a point not
  in original list, its skipped (as we want---no expansion of
  the original list). Earlier cases seem to still work with
  these changes. Time for the unit tests now. 

0.1.1
* Included a write-up from Rod that describes how the SNF is used to
  map points into "group coordinates". This is a key idea for the
  O(N) algorithm, which uses the mapping as a hash function for
  the kpoint list.

0.1.0
* Added functionality for an arbitrary shift. Seems to be working
  (tried several cases by hand). Next thing we need is a large
  collection of unit tests.
>>>>>>> ff2e3fa686db555283a0a33391bc85f6dd5bfbe3

<<< Note, still need to add functionality for skipping rotations that move a k-grid point
off the k-grid lattice. This can happen for strange shifts---which we should avoid, but it
would be nice for the code to do the right thing anyhow. >>> 

0.0.3
<<<<<<< HEAD

* Somehow, and old, incomplete copy got pushed and v. 0.0.2, if it was complete, got lost.
  This version seems to be working again, with the same functionality.
=======
* Somehow, and old, incomplete copy got pushed and v. 0.0.2, if it was
  complete, got lost. This version seems to be working again,
  with the same functionality.
>>>>>>> ff2e3fa686db555283a0a33391bc85f6dd5bfbe3

0.0.2
* Added files from laptop. First working copy.

0.0.1
<<<<<<< HEAD
* Initial import of rough, working copy of the folding code. Relies on routines in symlib
  (getPointGroup, SNF, HNF, etc.)

*Somehow, an old, incomplete copy got pushed and v. 0.0.2, if it was complete, got lost.
 This version seems to be working again, with the same functionality.
=======
* Initial import of rough, working copy of the folding code. Relies
  on routines in symlib (getPointGroup, SNF, HNF, etc.)
>>>>>>> ff2e3fa686db555283a0a33391bc85f6dd5bfbe3
