#pragma once
#include "FECore/FENLConstraint.h"
#include "FECore/FESolidDomain.h"
#include "FEPreStrainTransIsoMR.h"

//-----------------------------------------------------------------------------
class FEInSituStretch : public FENLConstraint
{
public:
	FEInSituStretch(FEModel* pfem);

public:
	bool Init();
	bool Augment(int naug);

	void Residual(FEGlobalVector& R) {};
	void StiffnessMatrix(FESolver* psolver) {};
	void Serialize(DumpFile& ar) {};
	virtual void Update() {}

private:
	bool CheckAugment(FESolidDomain* pdom, FEPreStrainTransIsoMR* pmat, int n);
	void DoAugment(FESolidDomain* pdom, FEPreStrainTransIsoMR* pmat, int n);

public:
	double	m_ltol;	//!< augmented Lagrangian tolerance

	// declare parameter list
	DECLARE_PARAMETER_LIST();
};
