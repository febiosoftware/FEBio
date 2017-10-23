#pragma once
#include <FECore/FENLConstraint.h>

//-----------------------------------------------------------------------------
// This class implements a constraint that enforces the distance between two nodes
class FEDistanceConstraint : public FENLConstraint
{
public:
	//! constructor
	FEDistanceConstraint(FEModel* pfem);

	bool Init();
	void Activate();
	void Residual(FEGlobalVector& R, const FETimeInfo& tp);
	void StiffnessMatrix(FESolver* psolver, const FETimeInfo& tp);
	bool Augment(int naug, const FETimeInfo& tp);
	void Serialize(DumpStream& ar);

	//! build connectivity for matrix profile
	void BuildMatrixProfile(FEGlobalMatrix& M);

	// update state
	void Reset();
	void Update(int niter, const FETimeInfo& tp);

	//! define matrix profile for constraint
	void BuildMatrixProfile();

public:
	double	m_eps;		//!< penalty parameter
	double	m_atol;		//!< augmented Lagrangian tolerance
	bool	m_blaugon;	//!< augmentation flag
	int		m_node[2];	//!< the two nodes that are connected
	int		m_nminaug;	//!< min number of augmentations
	int		m_nmaxaug;	//!< max number of augmentations

	double	m_l0;		//!< reference length
	double	m_Lm;		//!< Lagrange multiplier

	int	m_dofX;
	int	m_dofY;
	int	m_dofZ;

	DECLARE_PARAMETER_LIST();
};
