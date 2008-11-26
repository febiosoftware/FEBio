#pragma once

#include "FECore/vector.h"
#include "FECore/matrix.h"
#include <list>
using namespace std;

//-----------------------------------------------------------------------------
//! linear constraint enforced using augmented lagrangian

class FEAugLagLinearConstraint
{
public:
	class DOF
	{
	public:
		DOF() { node = bc = neq = -1; }
	public:
		int	node;	// the node to which this dof belongs to
		int	bc;		// the degree of freedom
		int	neq;	// the equation number (or -1 if none)
	};

	class SlaveDOF : public DOF
	{
	public:
		SlaveDOF() : val(0){}
		double	val;	// coefficient value
	};

public:
	FEAugLagLinearConstraint();

	void Residual(vector<double>& R);

	void Stiffness(vector<int>& en, vector<int>& elm, matrix& ke);

	bool Augment();

public:
	DOF				m_master;	// master degree of freedom
	list<SlaveDOF>	m_slave;	// list of slave nodes

	double	m_lam;	// lagrange multiplier
	double	m_eps;	// penalty parameter
	double	m_tol;	// convergence tolerance
};
