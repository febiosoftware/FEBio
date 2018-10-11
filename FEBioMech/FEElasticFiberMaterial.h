#pragma once
#include "FEElasticMaterial.h"
#include <FECore/FEVectorGenerator.h>

//-----------------------------------------------------------------------------
//! Base class for single fiber response

class FEElasticFiberMaterial : public FEElasticMaterial
{
public:
    FEElasticFiberMaterial(FEModel* pfem);

	bool Init() override;

	vec3d GetFiberVector(FEMaterialPoint& mp);

public:
	FEVectorGenerator*	m_fiberGenerator;

	DECLARE_FECORE_CLASS();
};
