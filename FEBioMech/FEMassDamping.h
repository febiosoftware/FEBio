#pragma once
#include <FEBioMech/FEBodyForce.h>

class FEMassDamping : public FEBodyForce
{
public:
	FEMassDamping(FEModel* fem);

	//! calculate the body force at a material point
	vec3d force(FEMaterialPoint& pt) override;

	//! calculate constribution to stiffness matrix
	mat3ds stiffness(FEMaterialPoint& pt) override;

private:
	double	m_C;

	DECLARE_FECORE_CLASS();
};
