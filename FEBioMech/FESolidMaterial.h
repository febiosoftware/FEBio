#pragma once
#include "FECore/FEMaterial.h"

//-----------------------------------------------------------------------------
//! Base class for solid-materials.
//! These materials need to define the stress and tangent functions.
//!
class FESolidMaterial : public FEMaterial
{
public:
	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt) = 0;

	//! calculate tangent stiffness at material point
	virtual tens4ds Tangent(FEMaterialPoint& pt) = 0;

	//! calculate derivative of stress w.r.t. solute concentration at material point
	mat3ds Tangent_Concentration(FEMaterialPoint& pt, const int isol);
	
	//! return the material density
	virtual double Density() { return m_density; }

	//! return the molar mass
	virtual double MolarMass() { return m_molarmass; }

	//! calculate strain energy density at material point \todo remove this
	virtual double StrainEnergy(FEMaterialPoint& pt) { return 0; }

	//! calculate tangent of strain energy density with solid density at material point \todo remove this
	virtual double Tangent_SE_Density(FEMaterialPoint& pt) { return 0;}

	//! calculate tangent of stress with solid density at material point \todo remove this
	virtual mat3ds Tangent_Stress_Density(FEMaterialPoint& pt) { return mat3ds(); }

protected:
	double	m_density;	//!< material density
	double	m_molarmass;//!< material molar mass (molecular weight)
};
