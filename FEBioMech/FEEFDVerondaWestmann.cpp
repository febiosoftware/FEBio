#include "stdafx.h"
#include "FEEFDVerondaWestmann.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEEFDVerondaWestmann, FEUncoupledMaterial)
	ADD_PARAMETER(m_VW.m_c1, "c1");
	ADD_PARAMETER(m_VW.m_c2, "c2");
	ADD_PARAMETER(m_EFD.m_beta, 3, "beta");
	ADD_PARAMETER(m_EFD.m_ksi , 3, "ksi" );
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEEFDVerondaWestmann::FEEFDVerondaWestmann(FEModel* pfem) : FEUncoupledMaterial(pfem), m_VW(pfem), m_EFD(pfem) 
{
	// although we don't use K of the child materials, we need to set it to a non-zero value
	// otherwise FEBio will complain
	m_VW.m_K = 1.0;
	m_EFD.m_K = 1.0;
}

//-----------------------------------------------------------------------------
bool FEEFDVerondaWestmann::Init()
{
	if (FEUncoupledMaterial::Init() == false) return false;
	if (m_VW.Init() == false) return false;
	if (m_EFD.Init() == false) return false;
	return true;
}

//-----------------------------------------------------------------------------
void FEEFDVerondaWestmann::Serialize(DumpStream& ar)
{
	FEUncoupledMaterial::Serialize(ar);
	m_VW.Serialize(ar);
	m_EFD.Serialize(ar);
}

//-----------------------------------------------------------------------------
mat3ds FEEFDVerondaWestmann::DevStress(FEMaterialPoint& pt)
{
	return m_VW.DevStress(pt) + m_EFD.DevStress(pt);
}

//-----------------------------------------------------------------------------
tens4ds FEEFDVerondaWestmann::DevTangent(FEMaterialPoint& pt)
{
	return m_VW.DevTangent(pt) + m_EFD.DevTangent(pt);
}

//-----------------------------------------------------------------------------
//! calculate deviatoric strain energy density
double FEEFDVerondaWestmann::DevStrainEnergyDensity(FEMaterialPoint& pt)
{
    return m_VW.DevStrainEnergyDensity(pt) + m_EFD.DevStrainEnergyDensity(pt);
}

