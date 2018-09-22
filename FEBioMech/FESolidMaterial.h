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
	FESolidMaterial(FEModel* pfem);

	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt) = 0;

	//! calculate tangent stiffness at material point
	virtual tens4ds Tangent(FEMaterialPoint& pt) = 0;

	//! calculate the 2nd Piola-Kirchhoff stress at material point
	virtual mat3ds PK2Stress(FEMaterialPoint& pt, const mat3ds E);

	//! calculate material tangent stiffness at material point
	virtual tens4ds MaterialTangent(FEMaterialPoint& pt, const mat3ds E);

	//! return the material density
	virtual double Density();

	//! return the material density
	virtual void SetDensity(const double d);

	//! Is this a rigid material or not
	virtual bool IsRigid() const { return false; }
    
protected:
	double	m_density;	//!< material density

	DECLARE_PARAMETER_LIST();
};
