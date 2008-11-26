#include "StdAfx.h"
#include "FEAugLagLinearConstraint.h"

FEAugLagLinearConstraint::FEAugLagLinearConstraint()
{
	m_lam = 0;
	m_eps = 1;
	m_tol = 0.1;
}

bool FEAugLagLinearConstraint::Augment()
{
	return true;
}

void FEAugLagLinearConstraint::Residual(vector<double> &R)
{
	
}

void FEAugLagLinearConstraint::Stiffness(vector<int>& en, vector<int>& elm, matrix& ke)
{

}
