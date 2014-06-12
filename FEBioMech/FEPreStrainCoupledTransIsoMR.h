#pragma once
#include "FEPreStrainTransIsoMR.h"
#include "FECoupledTransIsoMooneyRivlin.h"

//-----------------------------------------------------------------------------
class FEPreStrainCoupledTransIsoMR: public FEElasticMaterial
{
public:
	//! constructor
	FEPreStrainCoupledTransIsoMR(FEModel* pfem);

public:
	double	m_ltrg;	//!< target stretch

public:
	//! create material point data for this material
	FEMaterialPoint* CreateMaterialPointData();

	//! calculate deviatoric stress at material point
	mat3ds Stress(FEMaterialPoint& pt);

	//! calculate deviatoric tangent stiffness at material point
	tens4ds Tangent(FEMaterialPoint& pt);

	//! target fiber stretch
	double FiberStretch(FEMaterialPoint& pt);

protected:
	// calculate the prestrain deformation gradient
	mat3d PreStrainDeformationGradient(FEMaterialPoint& mp);

	// declare parameter list
	DECLARE_PARAMETER_LIST();

private:
	FECoupledTransIsoMooneyRivlin	m_mat;
};
