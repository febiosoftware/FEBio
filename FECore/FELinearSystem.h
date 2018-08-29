#pragma once
#include "FEGlobalMatrix.h"
#include "matrix.h"
#include <vector>
using namespace std;

//-----------------------------------------------------------------------------
// Experimental class to see if all the assembly operations can be moved to a class
// and out of the solver class.
class FECORE_API FELinearSystem
{
public:
	// Constructor
	// Takes a FEGlobalMatrix class K that will store the actual stiffness matrix
	// and a vector F which contains the assembled contribution of the prescribed 
	// degrees of freedom. The F vector must be added to the "force" vector. The u 
	// vector contains the nodal values of the prescribed degrees of freedom.
	FELinearSystem(FEGlobalMatrix& K, vector<double>& F, vector<double>& u);

public:
	// Assembly routine
	// This assembles the element stiffness matrix ke into the global matrix.
	// The contributions of prescribed degrees of freedom will be stored in m_F
	void AssembleLHS(vector<int>& lm, matrix& ke);

	// This assembles a matrix to the RHS by pre-multiplying the matrix with the 
	// prescribed value array U and then adding it to F
	void AssembleRHS(vector<int>& lm, matrix& ke, vector<double>& U);

	// This assembles a vetor to the RHS
	void AssembleRHS(vector<int>& lm, vector<double>& fe);

private:
	FEGlobalMatrix& m_K;
	vector<double>&	m_F;
	vector<double>&	m_u;
};
