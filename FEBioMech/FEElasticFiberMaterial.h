#pragma once
#include "FEElasticMaterial.h"

//-----------------------------------------------------------------------------
//! Base class for single fiber response

class FEElasticFiberMaterial : public FEElasticMaterial
{
public:
    FEElasticFiberMaterial(FEModel* pfem);

	FEMaterialPoint* CreateMaterialPointData() override;

	vec3d GetFiberVector(FEMaterialPoint& mp);

protected:
	// NOTE: Some fiber materials define a theta, phi parameter to define the fiber vector.
	//       Although this is deprecated, for backward compatibility this was feature was moved here
	double	m_thd;
	double	m_phd;

	DECLARE_FECORE_CLASS();
};
