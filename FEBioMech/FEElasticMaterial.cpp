#include "stdafx.h"
#include "FEElasticMaterial.h"
#include "FECore/FEModel.h"

BEGIN_FECORE_CLASS(FEElasticMaterial, FESolidMaterial)
END_FECORE_CLASS();

FEElasticMaterial::FEElasticMaterial(FEModel* pfem) : FESolidMaterial(pfem)
{ 
	m_density = 1;
	AddDomainParameter(new FEElasticStress());
}

//-----------------------------------------------------------------------------
FEElasticMaterial::~FEElasticMaterial()
{ 
	
}

//-----------------------------------------------------------------------------
//! return the strain energy density
double FEElasticMaterial::StrainEnergyDensity(FEMaterialPoint& pt) { return 0; }

//-----------------------------------------------------------------------------
FEElasticStress::FEElasticStress() : FEDomainParameter("stress")
{

}

FEParamValue FEElasticStress::value(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();
	return ep.m_s;
}
