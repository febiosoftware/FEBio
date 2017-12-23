#pragma once

#include "FECore/vector.h"
#include "FECore/matrix.h"
#include "FECore/FENLConstraint.h"
#include <list>
using namespace std;

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
		DOF() { node = bc = 0; val = 0.0; }
	public:
		int	node;		// the node to which this dof belongs to
		int	bc;			// the degree of freedom
		double	val;	// coefficient value
	};

	typedef list<FEAugLagLinearConstraint::DOF>::iterator Iterator;

public:
	//! constructor
	FEAugLagLinearConstraint() { m_lam = 0; }

	//! serialize data to archive
	void Serialize(DumpStream& ar);

public:
	list<DOF>	m_dof;	//!< list of participating dofs
	double		m_lam;	//!< lagrange multiplier
};

//-----------------------------------------------------------------------------
//! This class manages a group of linear constraints

class FELinearConstraintSet : public FENLConstraint
{
public:
	//! constructor
	FELinearConstraintSet(FEModel* pfem);

	//! add a linear constraint to the list
	void add(FEAugLagLinearConstraint* plc) { m_LC.push_back(plc); }

public:
	//! serialize data to archive
	void Serialize(DumpStream& ar) override;

	//! add the linear constraint contributions to the residual
	void Residual(FEGlobalVector& R, const FETimeInfo& tp) override;

	//! add the linear constraint contributions to the stiffness matrix
	void StiffnessMatrix(FESolver* psolver, const FETimeInfo& tp) override;

	//! do the augmentation
	bool Augment(int naug, const FETimeInfo& tp) override;

	//! build connectivity for matrix profile
	void BuildMatrixProfile(FEGlobalMatrix& M) override;

protected:
	//! calculate the constraint value
	double constraint(FEAugLagLinearConstraint& LC);

public:
	list<FEAugLagLinearConstraint*>	m_LC;	//!< list of linear constraints

public:
	double	m_tol;	//!< augmentation tolerance
	double	m_eps;	//!< penalty factor
    double  m_rhs;  //!< right-hand-side of linear constraint equation
	int		m_naugmax;	//!< max nr of augmentations
	int		m_naugmin;	//!< min nf of augmentations

	int	m_nID;		//!< ID of manager

	DECLARE_PARAMETER_LIST();
};
