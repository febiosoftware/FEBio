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
	void Residual(FEGlobalVector& R);
	void StiffnessMatrix(FESolver* psolver);
	bool Augment(int naug);
	void Serialize(DumpFile& ar);
	void ShallowCopy(DumpStream& dmp, bool bsave);


	// update state
	void Reset();
	void Update();

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

	DECLARE_PARAMETER_LIST();
};
