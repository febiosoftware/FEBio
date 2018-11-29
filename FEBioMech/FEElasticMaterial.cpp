#include "stdafx.h"
#include "FEElasticMaterial.h"
#include "FECore/FEModel.h"

BEGIN_FECORE_CLASS(FEElasticMaterial, FESolidMaterial)
END_FECORE_CLASS();

FEElasticMaterial::FEElasticMaterial(FEModel* pfem) : FESolidMaterial(pfem)
{ 
	m_density = 1; 
}

//-----------------------------------------------------------------------------
FEElasticMaterial::~FEElasticMaterial()
{ 
	
}

//-----------------------------------------------------------------------------
//! return the strain energy density
double FEElasticMaterial::StrainEnergyDensity(FEMaterialPoint& pt) { return 0; }
