#include "stdafx.h"
#include "FESolidMaterial.h"

// Material parameters for FEElasticMaterial
BEGIN_PARAMETER_LIST(FESolidMaterial, FEMaterial)
	ADD_PARAMETER(m_density, FE_PARAM_DOUBLE, "density");
END_PARAMETER_LIST();

FESolidMaterial::FESolidMaterial(FEModel* pfem) : FEMaterial(pfem) {}

//! return the material density
double FESolidMaterial::Density() { return m_density; }
