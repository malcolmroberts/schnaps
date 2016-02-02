SCHNAPS
=======

Solveur pour les lois de Conservation Hyperboliques Non-linéaires
Appliqué aux PlasmaS

Mode d'emploi

téléchargement (nécessite git)

Accès développeur:
git clone git+ssh://<gforge_account_name>\@scm.gforge.inria.fr//gitroot/schnaps/schnaps.git

Accès lecture seule:
git clone https://gforge.inria.fr/git/schnaps/schnaps.git

se placer dans le dossier SCHNAPS

puis:

cmake .

make

ctest

une série de tests démarre. Il est conseillé de faire repasser ces
tests après toute modification du code.

Pour lancer schnaps:

générer le maillage:

gmsh disque.geo -3

./schnaps

Cet exemple consiste à résoudre l'équation de transport à
l'intérieur d'un cylindre. La visualisation  des résultats
utilise gmsh

gmsh dgvisu.msh

(puis tools -> options -> view[0] -> adapt visualization grid pour
afficher une image plus jolie...)

Adapter le fichier source "schnaps.c" pour traiter des cas avec
d'autres maillages. Des exemples de maillages se trouvent dans le
dossier geo.

Génération de la documentation (avec doxygen):

cd doc/
doxygen doxyschnaps

SCHNAPS est sous licence CeCILL:

http://www.cecill.info/licences/Licence_CeCILL_V1.1-US.html

Génération de la documentation (avec doxygen):

cd doc/
doxygen doxyschnaps
 *
 *
 */

Installation de StarPU:

Télécharger FxT
mkdir /usr/local/fxtdir
cd FxT
./bootstrap
export FXTDIR=/usr/local/fxtdir
./configure --prefix=$FXTDIR
make -j4
make install

Télécharger StarPU a partir de la base svn (voir site StarPU)
brew 
cd StarPu/
mkdir build
./autogen.sh
(installer les paquets brew manquants)
../configure --with-fxt=$FXTDIR
make -j4
make check (facultatif)
make install

Désactiver les codelets opencl (éventuellement)
export STARPU_NOPENCL=0

Télécharger VITE a partir de la base svn (voir site vite trace)
brew install graphviz


