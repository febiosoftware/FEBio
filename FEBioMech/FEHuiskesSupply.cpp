#include "stdafx.h"
#include "FEHuiskesSupply.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_PARAMETER_LIST(FEHuiskesSupply, FESolidSupply)
ADD_PARAMETER(m_B, FE_PARAM_DOUBLE, "B");
ADD_PARAMETER(m_k, FE_PARAM_DOUBLE, "k");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Constructor. 
FEHuiskesSupply::FEHuiskesSupply(FEModel* pfem) : FESolidSupply(pfem)
{
	m_B = m_k = 0;
}

//-----------------------------------------------------------------------------
//! Initialization. 
void FEHuiskesSupply::Init()
{
}

//-----------------------------------------------------------------------------
//! Solid supply
double FEHuiskesSupply::Supply(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	FERemodelingMaterialPoint& rpt = *mp.ExtractData<FERemodelingMaterialPoint>();
	double rhor = pt.m_rhor;
	double sed = rpt.m_sed;
	double rhorhat = m_B*(sed/rhor - m_k);
	return rhorhat;
}

//-----------------------------------------------------------------------------
//! Tangent of solid supply with respect to strain
mat3ds FEHuiskesSupply::Tangent_Supply_Strain(FEMaterialPoint &mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    mat3ds ruhat = pt.m_s*(m_B/pt.m_rhor);
	return ruhat;
}

//-----------------------------------------------------------------------------
//! Tangent of solid supply with respect to referential density
double FEHuiskesSupply::Tangent_Supply_Density(FEMaterialPoint &mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	FERemodelingMaterialPoint& rpt = *mp.ExtractData<FERemodelingMaterialPoint>();
    double rhor = pt.m_rhor;
    double sed = rpt.m_sed;
    double dsed = rpt.dsed;
	return (dsed - sed/rhor)*m_B/rhor;
}

