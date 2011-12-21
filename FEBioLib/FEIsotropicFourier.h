#pragma once
#include "FEHeatTransferMaterial.h"

//-----------------------------------------------------------------------------
// Isotropic Fourer heat-transfer material
class FEIsotropicFourier : public FEHeatTransferMaterial
{
public:
	FEIsotropicFourier() {}
	void Init();

public:
	double	m_k;	//!< heat conductivity
	double	m_c;	//!< heat capacitance
	double	m_rho;	//!< density

public:
	void Conductivity(double D[3][3]);

	// declare parameter list
	DECLARE_PARAMETER_LIST();
};
