#include "stdafx.h"
#include "FEHuiskesSupply.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEHuiskesSupply, FESolidSupply)
	ADD_PARAMETER(m_B, "B");
	ADD_PARAMETER(m_k, "k");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor. 
FEHuiskesSupply::FEHuiskesSupply(FEModel* pfem) : FESolidSupply(pfem)
{
	m_B = m_k = 0;
}

//-----------------------------------------------------------------------------
//! Solid supply
double FEHuiskesSupply::Supply(FEMaterialPoint& mp)
{
	FERemodelingMaterialPoint& rpt = *mp.ExtractData<FERemodelingMaterialPoint>();
	double rhor = rpt.m_rhor;
	double sed = rpt.m_sed;
	double rhorhat = m_B*(sed/rhor - m_k);
	return rhorhat;
}

//-----------------------------------------------------------------------------
//! Tangent of solid supply with respect to strain
mat3ds FEHuiskesSupply::Tangent_Supply_Strain(FEMaterialPoint &mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	FERemodelingMaterialPoint& rpt = *mp.ExtractData<FERemodelingMaterialPoint>();
    mat3ds ruhat = pt.m_s*(m_B/rpt.m_rhor);
	return ruhat;
}

//-----------------------------------------------------------------------------
//! Tangent of solid supply with respect to referential density
double FEHuiskesSupply::Tangent_Supply_Density(FEMaterialPoint &mp)
{
	FERemodelingMaterialPoint& rpt = *mp.ExtractData<FERemodelingMaterialPoint>();
    double rhor = rpt.m_rhor;
    double sed = rpt.m_sed;
    double dsed = rpt.m_dsed;
	return (dsed - sed/rhor)*m_B/rhor;
}

