/* Sample program for the UMFPACK sparse matrix solver
*/

#include <stdio.h>
#include "solverumfpack.h"

int smalltestumfpack(void)
{
	int n = 5;
	double x[5];
	void *Symbolic, *Numeric;
	int i;

	/* cumulative count of entries, as matrix is scanned columnwise */
	int Ap[] = { 0, 2, 5, 9, 10, 12 };

	/* row indices of entries, as matrix is scanned columnwise */
	int Ai[] =    { 0, 1, 0,  2, 4, 1,  2, 3, 4, 2, 1, 4 };

	/* matrix entries */
	double Ax[] = { 2, 3, 3, -1, 4, 4, -3, 1, 2, 2, 6, 1 };

	/* the right hand side */
	double b[] = { 8, 45, -3, 3, 19 };

	/* symbolic analysis */
	umfpack_di_symbolic(n, n, Ap, Ai, Ax, &Symbolic, NULL, NULL);

	/* LU factorization */
	umfpack_di_numeric(Ap, Ai, Ax, Symbolic, &Numeric, NULL, NULL);
	umfpack_di_free_symbolic(&Symbolic);

	/* solve system */
	umfpack_di_solve(UMFPACK_A, Ap, Ai, Ax, x, b, Numeric, NULL, NULL);
	umfpack_di_free_numeric(&Numeric);

	for (i = 0; i < n; i++)
		printf("x[%d] = %g\n", i, x[i]);

	return 0;
}