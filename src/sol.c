/* sol.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

//#include "f2c.h"
#include <stdio.h>

/* Table of constant values */

static int c__1 = 1;

/* Subroutine */ int sol_(double *vkgs, double *vkgd, double *
	vkgi, double *vfg, int *kld, double *vu, int *neq, 
	int *mp, int *ifac, int *isol, int *nsym, double *
	energ, int *ier)
{
    /* Initialized data */

    static double vzero = 0.;

    /* Format strings */
    static char fmt_8000[] = "(\002 * sol pivot nul equation\002,i5)";

    /* System generated locals */
    int i__1, i__2, i__3, i__4;

    /* Builtin functions */
    //int s_wsfe(cilist *), do_fio(int *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static int i__;
    static double c1, c2;
    static int j1, j2, ic, ij, ik, jbk, jck, jhj, jhk, lhk, jhj1, jhk1, 
	    lhk1;
    extern double scal_(double *, double *, int *);
    static int imin, imax, imin1;
    static double cdiag;

    /* Fortran I/O blocks */
    //static cilist io___22 = { 0, 0, 0, fmt_8000, 0 };


/*   resolution d'un systeme lineaire symetrique ou non. la matrice est */
/*   stockee par ligne de ciel,en memoire dans les tables vkgs,vkgd,vkgi */

/*       entrees */
/*          vkgs,vkgd,vkgi    matrice du systeme : parties superieure, */
/*                            diagonale, inferieure (double precision) */
/*          vfg               second membre */
/*          kld               pointeurs vers les hauts de colonne */
/*          vu                vecteur solution (qui peut etre vfg) */
/*          neq               nombre d'equations */
/*          mp                unite logique d'impression */
/*          ifac              si ifac.eq.1 triangularisation de */
/*                            la matrice */
/*          isol              si isol.eq.1 calcul de la solution a */
/*                            partir de la matrice triangularisee */
/*          nsym              indice de probleme non symetrique */
/*       sorties */
/*          vkgs,vkgd,vkgi    matrice triangularisee (si ifac.eq.1) */
/*          vfg               solution (si isol.eq.1) */
/*          energ             energie du systeme (si nsym.eq.0) */
/*          ier               mis a 1 si pivot nul rencontre */

/* =========================== debut des declarations ==================== */
    /* Parameter adjustments */
    --vu;
    --kld;
    --vfg;
    --vkgi;
    --vkgd;
    --vkgs;

    /* Function Body */
/* =========================== debut du code executable ================== */

/* -------  traitement */

    ik = 1;
    if (vkgd[1] == vzero) {
	goto L800;
    }
    *energ = vzero;
    *ier = 0;
    if (*isol == 1) {
	i__1 = *neq;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    vu[i__] = vfg[i__];
	}
    }

/* -------  pour chaque colonne ik a modifier */

    jhk = 1;
    i__1 = *neq;
    for (ik = 2; ik <= i__1; ++ik) {

/* -------  pointeur du haut de la colonne suivante ik+1 */

	jhk1 = kld[ik + 1];

/* -------  hauteur de la colonne ik (hors termes superieur et diagonal) */

	lhk = jhk1 - jhk;
	lhk1 = lhk - 1;

/* -------  ligne du premier terme a modifier dans la colonne ik */

	imin = ik - lhk1;
	imin1 = imin - 1;

/* -------  ligne du dernier terme a modifier dans la colonne ik */

	imax = ik - 1;
	if (lhk1 < 0) {
	    goto L100;
	}
	if (*ifac != 1) {
	    goto L90;
	}
	if (*nsym == 1) {
	    vkgi[jhk] /= vkgd[imin1];
	}
	if (lhk1 == 0) {
	    goto L40;
	}

/* -------  modifier les termes non diagonaux de la colonne ik */

	jck = jhk + 1;
	jhj = kld[imin];

/* -------  pour chaque terme place en jck, correspondant a la colonne ij */

	i__2 = imax;
	for (ij = imin; ij <= i__2; ++ij) {
	    jhj1 = kld[ij + 1];

/* -------  nombre de termes modificatifs du terme place en jck */

/* Computing MIN */
	    i__3 = jck - jhk, i__4 = jhj1 - jhj;
	    ic = min(i__3,i__4);
	    if (ic <= 0 && *nsym == 0) {
		goto L20;
	    }
	    c1 = vzero;
	    if (ic <= 0) {
		goto L17;
	    }
	    j1 = jhj1 - ic;
	    j2 = jck - ic;
	    if (*nsym == 1) {
		goto L15;
	    }
	    vkgs[jck] -= scal_(&vkgs[j1], &vkgs[j2], &ic);
	    goto L20;
L15:
	    vkgs[jck] -= scal_(&vkgi[j1], &vkgs[j2], &ic);
	    c1 = scal_(&vkgs[j1], &vkgi[j2], &ic);
L17:
	    vkgi[jck] = (vkgi[jck] - c1) / vkgd[ij];
L20:
	    ++jck;
/* L30: */
	    jhj = jhj1;
	}

/* -------  modifier le terme diagonal */

L40:
	jck = jhk;
	cdiag = vzero;
	i__2 = imax;
	for (ij = imin1; ij <= i__2; ++ij) {
	    c1 = vkgs[jck];
	    if (*nsym == 1) {
		goto L50;
	    }
	    c2 = c1 / vkgd[ij];
	    vkgs[jck] = c2;
	    goto L60;
L50:
	    c2 = vkgi[jck];
L60:
	    cdiag += c1 * c2;
/* L70: */
	    ++jck;
	}
	vkgd[ik] -= cdiag;
	if (vkgd[ik] == 0.f) {
	    goto L800;
	}

/* -------  resolution du systeme triangulaire inferieur */

L90:
	if (*isol != 1) {
	    goto L100;
	}
	if (*nsym != 1) {
	    vu[ik] = vfg[ik] - scal_(&vkgs[jhk], &vu[imin1], &lhk);
	}
	if (*nsym == 1) {
	    vu[ik] = vfg[ik] - scal_(&vkgi[jhk], &vu[imin1], &lhk);
	}
L100:
	jhk = jhk1;
    }
    if (*isol != 1) {
	goto L9999;
    }

/* -------  resolution du systeme diagonal : */

    if (*nsym == 1) {
	goto L120;
    }
    i__1 = *neq;
    for (ik = 1; ik <= i__1; ++ik) {
	c1 = vkgd[ik];
	if (c1 == vzero) {
	    goto L800;
	}
	c2 = vu[ik] / c1;
	vu[ik] = c2;
/* L110: */
	*energ += c1 * c2 * c2;
    }

/* -------  resolution du systeme triangulaire superieur */

L120:
    ik = *neq + 1;
    jhk1 = kld[ik];
L130:
    --ik;
    if (*nsym == 1) {
	vu[ik] /= vkgd[ik];
    }
    if (ik == 1) {
	goto L9999;
    }
    c1 = vu[ik];
    jhk = kld[ik];
    jbk = jhk1 - 1;
    if (jhk > jbk) {
	goto L150;
    }
    ij = ik - jbk + jhk - 1;
    i__1 = jbk;
    for (jck = jhk; jck <= i__1; ++jck) {
	vu[ij] -= vkgs[jck] * c1;
/* L140: */
	++ij;
    }
L150:
    jhk1 = jhk;
    goto L130;

/* -------  erreurs */

L800:
    if (*mp != 0) {
	/* io___22.ciunit = *mp; */
      printf("%s\n",fmt_8000);
	/* s_wsfe(&io___22); */
	/* do_fio(&c__1, (char *)&ik, (ftnlen)sizeof(int)); */
	/* e_wsfe(); */
	
    }
    *ier = 1;
    goto L9999;

/* -------  fin */

L9999:
    return 0;
/* ===========================   fin du module sol    ================== */
} /* sol_ */

double scal_(double *x, double *y, int *n)
{
    /* Initialized data */

    static double zero = 0.;

    /* System generated locals */
    int i__1;
    double ret_val;

    /* Local variables */
    static int i__;

/* ======================================================================= */
/* calcul du produit scalaire */
/* ======================================================================= */
    /* Parameter adjustments */
    --y;
    --x;

    /* Function Body */
/* ----------------------------------------------------------------------- */
    ret_val = zero;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ret_val += x[i__] * y[i__];
    }
    return ret_val;
} /* scal_ */

