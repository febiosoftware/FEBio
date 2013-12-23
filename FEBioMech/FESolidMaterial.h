#pragma once
#include "FECore/FEMaterial.h"

//-----------------------------------------------------------------------------
//! Base class for solid-materials.
//! These materials need to define the stress and tangent functions.
//!
class FESolidMaterial : public FEMaterial
{
public:
	//! constructor
	FESolidMaterial(FEModel* pfem) : FEMaterial(pfem) {}

	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt) = 0;

	//! calculate tangent stiffness at material point
	virtual tens4ds Tangent(FEMaterialPoint& pt) = 0;

	//! calculate derivative of stress w.r.t. solute concentration at material point
	mat3ds Tangent_Concentration(FEMaterialPoint& pt, const int isol);
	
	//! return the material density
	virtual double Density() { return m_density; }

protected:
	double	m_density;	//!< material density
};
