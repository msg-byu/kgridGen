0.2.4
*-* Deleted a portion of the code in symmetryReduceKpointList that
&nbsp; was a copy of code in findKptIndex. This fixed a bug.
*-* Started enforcing new lines when line length exceeded 70
&nbsp; characters.
*-* Tidied up HISTORY.md and findKptIndex.
*-* Made it so the offset in generateFullKpointList always lies
&nbsp; within the first unit cell.
*-* Added a few more unit tests for simple cubic lattices.

0.2.3 
*-* Added new unit tests for body-centered, cubic lattices.

0.2.2
*-* Fixed finite precision issues related to eps.
*-* Added a few unit tests for simple cubic lattices.
*-* Added an additional fail safe test.

0.2.1 Bug fix: fixed finite precision error when checking if a rotated
&nbsp; kpoint is off the k-grid. Was using an epsilon check without
&nbsp; using abs() and when the values were negative, the check
&nbsp; incorrectly failed. Replaced the condition with the "equal"
&nbsp; function from numerical_utilities module. Should have been
&nbsp; using that all along anyway---I coded it up so that I wouldn't
&nbsp; make these kinds of mistakes!
&nbsp; Replaced the line:
&nbsp; ```if (any((matmul(invK,roKpt) - nint(matmul(invK,roKpt))) > eps)) then```
&nbsp; with
&nbsp; ```if (.not. equal(matmul(invK,roKpt), nint(matmul(invK,roKpt)), eps)) then```
      
0.2.0
*-* Bug fixes, redefinition of shift. Arbitrary shifts seem to do
&nbsp; the right thing. If a rotation maps a kpoint onto a point not
&nbsp; in original list, its skipped (as we want---no expansion of
&nbsp; the original list). Earlier cases seem to still work with
&nbsp; these changes. Time for the unit tests now. 

0.1.1
*-* Included a write-up from Rod that describes how the SNF is used to
&nbsp; map points into "group coordinates". This is a key idea for the
&nbsp; O(N) algorithm, which uses the mapping as a hash function for
&nbsp; the kpoint list.

0.1.0
*-* Added functionality for an arbitrary shift. Seems to be working
&nbsp; (tried several cases by hand). Next thing we need is a large
&nbsp; collection of unit tests.

<<< Note, still need to add functionality for skipping rotations that
move a k-grid point off the k-grid lattice. This can happen for
strange shifts---which we should avoid, but it would be nice for the
code to do the right thing anyhow. >>> 

0.0.3
*-*Somehow, and old, incomplete copy got pushed and v. 0.0.2, if it was
&nbsp; complete, got lost. This version seems to be working again,
&nbsp; with the same functionality.

0.0.2
*-*Added files from laptop. First working copy.

0.0.1
*-* Initial import of rough, working copy of the folding code. Relies
&nbsp; on routines in symlib (getPointGroup, SNF, HNF, etc.)
