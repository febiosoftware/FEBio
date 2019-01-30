#pragma once
#include <FECore/CompactUnSymmMatrix.h>

namespace NumCore
{
	// write a matrix to file
	bool write_hb(CompactMatrix& K, const char* szfile);

	// write a vector to file
	bool write_vector(const vector<double>& a, const char* szfile);

	// calculate inf-norm of inverse matrix (only works with CRSSparsMatrix(1))
	double inverse_infnorm(CompactMatrix* A);

	// calculate condition number of a CRSSparseMatrix(1)
	double conditionNumber(CRSSparseMatrix* A);

	// create a random vector
	void randomVector(vector<double>& R, double vmin = 0.0, double vmax = 1.0);

	// inf-norm of a vector
	double infNorm(const std::vector<double>& x);

} // namespace NumCore
