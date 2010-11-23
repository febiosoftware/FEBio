#pragma once

#include "FECore/vector.h"
#include "FECore/matrix.h"
#include "FEStiffnessMatrix.h"
#include "DumpFile.h"
#include <list>
using namespace std;

class FEM;

//-----------------------------------------------------------------------------
//! linear constraint enforced using augmented lagrangian

class FEAugLagLinearConstraint
{
public:
	// this class describes a degree of freedom (dof) that
	// participates in the linear constraint
	class DOF
	{
	public:
		DOF() { node = bc = neq = -1; val = 0; }
	public:
		int	node;		// the node to which this dof belongs to
		int	bc;			// the degree of freedom
		int	neq;		// the equation number (or -1 if none)
		double	val;	// coefficient value
	};

	typedef list<FEAugLagLinearConstraint::DOF>::iterator Iterator;

public:
	//! constructor
	FEAugLagLinearConstraint() { m_lam = 0; }

	//! serialize data to archive
	void Serialize(DumpFile& ar);

public:
	list<DOF>	m_dof;	//!< list of participating dofs
	double		m_lam;	//!< lagrange multiplier
};

//-----------------------------------------------------------------------------
//! This class manages a group of linear constraints

class FELinearConstraintSet
{
public:
	//! constructor
	FELinearConstraintSet(FEM* pfem);

	//! add a linear constraint to the list
	void add(FEAugLagLinearConstraint* plc) { m_LC.push_back(plc); }

	//! add the linear constraint contributions to the residual
	void Residual(vector<double>& R);

	//! add the linear constraint contributions to the stiffness matrix
	void Stiffness();

	//! do the augmentation
	bool Augment(int naug);

	//! initialization
	void Init();

	//! serialize data to archive
	void Serialize(DumpFile& ar);

protected:
	//! calculate the constraint value
	double constraint(FEAugLagLinearConstraint& LC);

public:
	list<FEAugLagLinearConstraint*>	m_LC;	//!< list of linear constraints

public:
	FEM*	m_pfem;	//!< pointer to FEM data

	double	m_tol;	//!< augmentation tolerance
	double	m_eps;	//!< penalty factor

	int	m_naugmax;	//!< max nr of augmenations

	int	m_nID;		//!< ID of manager
};
