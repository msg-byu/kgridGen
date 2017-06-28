0.1.1 Included a write-up from Rod that describes how the SNF is used to map points into "group
coordinates". This is a key idea for the O(N) algorithm, which using the mapping as a hash
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
