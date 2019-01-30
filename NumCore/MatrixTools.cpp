#include <FECore/CompactMatrix.h>
#include "MatrixTools.h"
#include "PardisoSolver.h"
#include <stdlib.h>

bool NumCore::write_hb(CompactMatrix& K, const char* szfile)
{
	// This currently assumes a CRS matrix
	FILE* fp = fopen(szfile, "wb");
	if (fp == 0) return false;

	int	symmFlag = K.isSymmetric();
	int offset = K.Offset();
	int rowFlag = K.isRowBased();
	int nrow = K.Rows();
	int ncol = K.Columns();
	int nsize = nrow;	// this is the reason we currently assume CRS
	int nnz = K.NonZeroes();
	fwrite(&symmFlag, sizeof(symmFlag), 1, fp);
	fwrite(&offset, sizeof(offset), 1, fp);
	fwrite(&rowFlag, sizeof(rowFlag), 1, fp);
	fwrite(&nrow, sizeof(nrow), 1, fp);
	fwrite(&ncol, sizeof(ncol), 1, fp);
	fwrite(&nnz, sizeof(nnz), 1, fp);
	fwrite(K.Pointers(), sizeof(int), nsize + 1, fp);
	fwrite(K.Indices(), sizeof(int), nnz, fp);
	fwrite(K.Values(), sizeof(double), nnz, fp);

	fclose(fp);

	return true;
}

// write a vector to file
bool NumCore::write_vector(const vector<double>& a, const char* szfile)
{
	FILE* fp = fopen(szfile, "wb");
	if (fp == 0) return false;

	int N = (int)a.size();
	fwrite(&N, sizeof(int), 1, fp);
	fwrite(&a[0], sizeof(double), N, fp);

	fclose(fp);

	return true;
}

// calculate inf-norm of inverse matrix (only works with CRSSparsMatrix(1))
double NumCore::inverse_infnorm(CompactMatrix* A)
{
	PardisoSolver solver(0);
	if (solver.SetSparseMatrix(A) == false) return 0.0;
	if (solver.PreProcess() == false) return 0.0;
	if (solver.Factor() == false) return 0.0;

	int N = A->Rows();
	vector<double> e(N, 0.0), x(N, 0.0);
	vector<double> s(N, 0.0);
	for (int i = 0; i < N; ++i)
	{
		// get the i-th column of the inverse matrix
		e[i] = 1.0;
		solver.BackSolve(&x[0], &e[0]);

		// add to net row sums
		for (int j = 0; j < N; ++j) s[j] += fabs(x[j]);

		// reset e
		e[i] = 0.0;

		if (i % 100 == 0)
			fprintf(stderr, "%.2lg%%\r", 100.0 *i / N);
	}

	// get the max row sum
	double smax = s[0];
	for (int i = 1; i < N; ++i) if (s[i] > smax) smax = s[i];

	return smax;
}

// calculate condition number of a CRSSparseMatrix(1)
double NumCore::conditionNumber(CRSSparseMatrix* A)
{
	double norm = A->infNorm();
	double inorm = inverse_infnorm(A);
	double c = norm*inorm;
	return c;
}

inline double frand() { return rand() / (double)RAND_MAX; }

void NumCore::randomVector(vector<double>& R, double vmin, double vmax)
{
	int neq = (int)R.size();
	for (int i = 0; i<neq; ++i)
	{
		R[i] = vmin + frand()*(vmax - vmin);
	}
}

// norm of a vector
double NumCore::infNorm(const std::vector<double>& x)
{
	double m = 0.0;
	for (size_t i = 0; i < x.size(); ++i)
	{
		double xi = fabs(x[i]);
		if (xi > m) m = xi;
	}
	return m;
}
