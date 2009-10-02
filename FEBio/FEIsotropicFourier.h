#pragma once
#include "FEMaterial.h"

class FEIsotropicFourier : public FEMaterial
{
public:
	FEIsotropicFourier() {}
	void Init();

public:
	double	m_k;	//!< heat conductivity
	double	m_c;	//!< heat capacitance
	double	m_rho;	//!< density

public:
	FEMaterialPoint* CreateMaterialPointData() { return new FEHeatMaterialPoint; }

public:
	void Conductivity(double D[3][3]);

	// declare as registered
	DECLARE_REGISTERED(FEIsotropicFourier);

	// declare parameter list
	DECLARE_PARAMETER_LIST();
};
