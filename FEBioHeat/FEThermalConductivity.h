#pragma once
#include <FECore/FEMaterial.h>

//-----------------------------------------------------------------------------
//! This class implements the thermal conductivity material property.
class FEThermalConductivity : public FEMaterial
{
public:
	FEThermalConductivity(FEModel* pfem) : FEMaterial(pfem) {}

	// evaluate conductivity
	virtual mat3ds Conductivity(FEMaterialPoint& mp) = 0;

	// evaluate the derivative of the conductivity with respect to C
	virtual tens4ds Tangent_Conductivity_Strain(FEMaterialPoint& mp) = 0;

	// evaluate the derivative of the conductivity with respect to temperature
	virtual mat3ds Tangent_Conductivity_Temperature(FEMaterialPoint& mp) = 0;
};

//-----------------------------------------------------------------------------
//! This class implements a conductivity tensor which is contant in the reference configuration
class FEConstReferenceThermalConductivity : public FEThermalConductivity
{
public:
	FEConstReferenceThermalConductivity(FEModel* pfem);

public:
	// evaluate conductivity
	mat3ds Conductivity(FEMaterialPoint& mp);

	// evaluate the derivative of the conductivity with respect to C
	tens4ds Tangent_Conductivity_Strain(FEMaterialPoint& mp);

	// evaluate the derivative of the conductivity with respect to temperature
	mat3ds Tangent_Conductivity_Temperature(FEMaterialPoint& mp);

public:
	double	m_k0;	//!< reference thermal conductivity
	double	m_wt;	//!< thermal relaxation parameter

	DECLARE_PARAMETER_LIST();
};

//-----------------------------------------------------------------------------
//! This class implements a conductivity tensor which is contant in the current configuration
class FEConstThermalConductivity : public FEThermalConductivity
{
public:
	FEConstThermalConductivity(FEModel* pfem);

	// initialization
	void Init();

public:
	// evaluate conductivity
	mat3ds Conductivity(FEMaterialPoint& mp);

	// evaluate the derivative of the conductivity with respect to C
	tens4ds Tangent_Conductivity_Strain(FEMaterialPoint& mp);

	// evaluate the derivative of the conductivity with respect to temperature
	mat3ds Tangent_Conductivity_Temperature(FEMaterialPoint& mp);

public:
	double	m_k0;	//!< reference thermal conductivity
	double	m_wt;	//!< thermal relaxation parameter

	DECLARE_PARAMETER_LIST();
};
