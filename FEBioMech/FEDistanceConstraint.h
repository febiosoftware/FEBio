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

	double	m_L0;		//!< reference length
	vec3d	m_Fc;		//!< contact force

	DECLARE_PARAMETER_LIST();
};
