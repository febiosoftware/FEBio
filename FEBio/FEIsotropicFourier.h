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
	double BulkModulus() { return 0; }
	double Density() { return 0; }

	mat3ds Stress(FEMaterialPoint& pt) { return mat3ds(); }
	tens4ds Tangent(FEMaterialPoint& pt) { return tens4ds(); }

public:
	void Conductivity(double D[3][3]);

	// declare as registered
	DECLARE_REGISTERED(FEIsotropicFourier);

	// declare parameter list
	DECLARE_PARAMETER_LIST();
};
