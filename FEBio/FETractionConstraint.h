#pragma once

#include "FECore/vector.h"
#include "FECore/matrix.h"
#include "FEStiffnessMatrix.h"
#include "Archive.h"
#include <list>
using namespace std;

class FEM;

//-----------------------------------------------------------------------------
//! traction constraint enforced using augmented lagrangian

class FETractionConstraint
{
public:
	// this class describes a degree of freedom (dof) that
	// participates in the traction constraint
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

	typedef list<FETractionConstraint::DOF>::iterator Iterator;

public:
	//! constructor
	FETractionConstraint() { m_lam = 0; }

	//! serialize data to archive
	void Serialize(Archive& ar);

public:
	list<DOF>	m_dof;	//!< list of participating dofs
	double		m_lam;	//!< lagrange multiplier
};

//-----------------------------------------------------------------------------
//! This class manages a group of traction constraints

class FETractionConstraintSet
{
public:
	//! constructor
	FETractionConstraintSet(FEM* pfem);

	//! add a traction constraint to the list
	void add(FETractionConstraint* prc) { m_RC.push_back(prc); }

	//! add the traction constraint contributions to the residual
	void Residual(vector<double>& R);

	//! add the traction constraint contributions to the stiffness matrix
	void Stiffness();

	//! do the augmentation
	bool Augment(int naug);

	//! initialization
	void Init();

	//! serialize data to archive
	void Serialize(Archive& ar);

protected:
	//! calculate the constraint value
	double constraint(FETractionConstraint& RC);

public:
	list<FETractionConstraint*>	m_RC;	//!< list of linear constraints

public:
	FEM*	m_pfem;	//!< pointer to FEM data

	double	m_tol;	//!< augmentation tolerance
	double	m_eps;	//!< penalty factor

	int	m_naugmax;	//!< max nr of augmenations

	int	m_nID;		//!< ID of manager
};
