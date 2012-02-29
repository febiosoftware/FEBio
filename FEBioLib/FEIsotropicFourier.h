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
	//! get the material's conductivity tensor
	void Conductivity(double D[3][3]);

	//! get the material's capacitance
	double Capacitance() { return m_c; }

	//! get the material's density
	double Density() { return m_rho; }

	// declare parameter list
	DECLARE_PARAMETER_LIST();
};
