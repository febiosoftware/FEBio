#pragma once
#include "FECore/FENLConstraint.h"

//-----------------------------------------------------------------------------
class FEInSituStretch : public FENLConstraint
{
public:
	FEInSituStretch(FEModel* pfem);

public:
	void Init();
	bool Augment(int naug);

	void Residual(FEGlobalVector& R) {};
	void StiffnessMatrix(FESolver* psolver) {};
	void Serialize(DumpFile& ar) {};
	virtual void Update() {}

public:
	double	m_ltol;	//!< augmented Lagrangian tolerance

	// declare parameter list
	DECLARE_PARAMETER_LIST();
};
