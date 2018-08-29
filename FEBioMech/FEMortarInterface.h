#pragma once
#include "FEContactInterface.h"
#include "FEMortarContactSurface.h"

//-----------------------------------------------------------------------------
// Base class for mortar-type contact formulations
class FEMortarInterface : public FEContactInterface
{
public:
	//! constructor
	FEMortarInterface(FEModel* pfem);

	//! update the mortar weights
	void UpdateMortarWeights(FESurface& ss, FESurface& ms);

	//! update the nodal gaps
	void UpdateNodalGaps(FEMortarContactSurface& ss, FEMortarContactSurface& ms);

protected:
	matrix	m_n1;	//!< integration weights n1_AB
	matrix	m_n2;	//!< integration weights n2_AB

private:
	// integration rule
	FESurfaceElementTraits*	m_pT;
};
