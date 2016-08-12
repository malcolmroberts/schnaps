SCHNAPS
=======

Solver for Conservative Hypebolic Non-linear systems Applied to PlasmaS

Solveur pour les lois de Conservation Hyperboliques Non-linéaires
Appliqué aux PlasmaS

Mode d'emploi

Downloads:

Develloper access:
git clone git+ssh://<gforge_account_name>\@scm.gforge.inria.fr//gitroot/schnaps/schnaps.git

Read-only access:
git clone https://gforge.inria.fr/git/schnaps/schnaps.git

se placer dans le dossier SCHNAPS

Compilation:

From the schnaps folder, first

mkdir build

cd build

cmake ..

make

gmsh ../geo/disque.geo -3

ctest

By default, the OpenCL code runs on platform 0, device 0.  This
can be changed in schnaps via command-line argument.  The unit tests
run on the default platform/device, which can be changed via

cmake -D_CL_PLATFORM=1 -D_CL_DEVICE=1  ..

By default, computations are done in double-precision.  This can be
changed via

cmake -DSINGLE_PRECISION:BOOL=ON

./schnaps

The main schnaps program is schnaps.  To generate .msh files from .geo
files, one run

gmsh <file>.geo -3

The resulting .msh file can be passed to schnaps via command-line
arguments (see ./schnaps -h for a list of options).

The main executable schnaps will output to a gmsh file "dgvisu.msh" if
the one calls it with ./schnaps -w 1.  By default no output is
written.

The default arguments for schnaps results in a simulation of 2D
transport in a disc.  The results can be viewed via

gmsh dgvisu.msh

The default gmsh visualizer is quite coarse, but the visualization can
be improved by the setting "Adapt visuzliation grid" in
tools -> options -> View [0].

Adapter le fichier source "schnaps.c" pour traiter des cas avec
d'autres maillages. Des exemples de maillages se trouvent dans le
dossier geo.

There is some basic documentation with doxygen, which can be generated via

cd doc/
doxygen doxyschnaps

SCHNAPS est sous licence CeCILL:

Debugging:

valgrind doesn't play terribly well with the AMD drivers.  This can be
ameliorated by using the provided suppression file, libamdocl.supp, via

valgrind --suppressions=libamdocl.supp <executable>


SCHNAPS is under the CeCILL license:
http://www.cecill.info/licences/Licence_CeCILL_V1.1-US.html
