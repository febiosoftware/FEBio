#pragma once
#include "FEGlobalMatrix.h"
#include "matrix.h"
#include <vector>
using namespace std;

//-----------------------------------------------------------------------------
// Experimental class to see if all the assembly operations can be moved to a class
// and out of the solver class.
class FEStiffnessMatrix
{
public:
	// Constructor
	// Takes a FEGlobalMatrix class K that will store the actual stiffness matrix
	// and a vector F which contains the assembled contribution of the prescribed 
	// degrees of freedom. The F vector must be added to the "force" vector. The u 
	// vector contains the nodal values of the prescribed degrees of freedom.
	FEStiffnessMatrix(FEGlobalMatrix& K, vector<double>& F, vector<double>& u);

public:
	// Assembly routine
	// This assembles the element stiffness matrix ke into the global matrix.
	// The contributions of prescribed degrees of freedom will be store in m_F
	void Assemble(vector<int>& en, vector<int>& lm, matrix& ke);

private:
	FEGlobalMatrix& m_K;
	vector<double>&	m_F;
	vector<double>&	m_u;
};
