#include "stdafx.h"
#include "FEStiffnessMatrix.h"

//-----------------------------------------------------------------------------
FEStiffnessMatrix::FEStiffnessMatrix(FEGlobalMatrix& K, vector<double>& F, vector<double>& u) : m_K(K), m_F(F), m_u(u)
{
}

//-----------------------------------------------------------------------------
//! assemble global stiffness matrix
void FEStiffnessMatrix::Assemble(vector<int>& en, vector<int>& elm, matrix& ke)
{
	// assemble into the global stiffness
	m_K.Assemble(ke, elm);

	// check the prescribed contributions
	SparseMatrix& K = m_K;
	int N = ke.rows();
	int neq = m_K.Rows();

	// loop over columns
	for (int j = 0; j<N; ++j)
	{
		int J = -elm[j] - 2;
		if ((J >= 0) && (J<neq))
		{
			// dof j is a prescribed degree of freedom

			// loop over rows
			for (int i = 0; i<N; ++i)
			{
				int I = elm[i];
				if (I >= 0)
				{
					// dof i is not a prescribed degree of freedom
					m_F[I] -= ke[i][j] * m_u[J];
				}
			}

			// set the diagonal element of K to 1
			K.set(J, J, 1);

			// set the rhs vector to the prescribed value
			// that way the solution vector will contain the prescribed value
			m_F[J] = m_u[J];
		}
	}
}
