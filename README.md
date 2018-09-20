# kgridGen
Generator for k-point grids.

This code takes a general, integer multiple (3x3 matrix) of the reciprocal lattice vectors and
generates a corresponding integration grid. The code can also generate the symmetrically distinct
kpoints of the full list, with their corresponding weights. Arbitrary shifts of the k-grid
generating vectors are also allowed.

## To compile the FORTRAN library.

First we need a clean copy of `symlib` (https://github.com/msg-byu/symlib/).

```
git clone git@github.com:msg-byu/symlib.git
cd symlib/src/
make F90=gfortran
cd ../../
```

You can alternatively compile with the ifort compiler. New we can
compile the src code for `kgridGen`.

```
git clone git@github.com:msg-byu/kgridGen.git
cd kgridGen/src/
make F90=gfortran kpoints.x
```

This will compile all the needed libraries and the executable for the
kpoints code. Please note that the input for the code is in the driver
so in order to run different cases.

## To compile the python libraries.

To compile the python wrapped version first complete the steps needed
for compiling the fortran code. Then do the following in the
`kgridGen/wrap` directory.

```
pip install f90wrap
make F90=gfortran
```

This will create the shared libraries `_kpoints.so` and `_kgridgen.so`
which will need to be added to the python path or kept in the folder
with the kpoints.py file.


For MSG students: 
Bounty for bugs: $20
Bounty for fixes: $20
