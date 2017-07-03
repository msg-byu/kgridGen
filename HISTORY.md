0.2.3 Added new unit tests for body-centered, cubic lattices.

0.2.2 Fixed finite precision issues related to eps. Added a few unit tests for simple cubic
      lattices. Added an additional fail safe test.
      
0.2.1 Bug fix: fixed finite precision error when checking if a rotated kpoint is off the k-grid. Was
using an epsilon check without using abs() and when the values were negative, the check incorrectly
failed. Replaced the condition with the "equal" function from numerical_utilities module. Should
have been using that all along anyway---I coded it up so that I wouldn't make these kinds of
mistakes!

The line:
        if (any((matmul(invK,roKpt) - nint(matmul(invK,roKpt))) > eps)) then
with
       if (.not. equal(matmul(invK,roKpt), nint(matmul(invK,roKpt)), eps)) then

0.2.0 Bug fixes, redefinition of shift. Arbitrary shifts seem to do the right thing. If a rotation
maps a kpoint onto a point not in original list, its skipped (as we want---no expansion of the
original list). Earlier cases seem to still work with these changes. Time for the unit tests now. 

0.1.1 Included a write-up from Rod that describes how the SNF is used to map points into "group
coordinates". This is a key idea for the O(N) algorithm, which uses the mapping as a hash
function for the kpoint list.

0.1.0 Added functionality for an arbitrary shift. Seems to be working (tried several cases by
hand). Next thing we need is a large collection of unit tests.

<<< Note, still need to add functionality for skipping rotations that move a k-grid point off the
k-grid lattice. This can happen for strange shifts---which we should avoid, but it would be nice for
the code to do the right thing anyhow. >>> 

0.0.3 Somehow, and old, incomplete copy got pushed and v. 0.0.2, if it was complete, got lost. This
version seems to be working again, with the same functionality.

0.0.2 Added files from laptop. First working copy.

0.0.1 Initial import of rough, working copy of the folding code. Relies on routines in symlib
(getPointGroup, SNF, HNF, etc.)
