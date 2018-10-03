#include "stdafx.h"
#include "FEEFDMooneyRivlin.h"

// define the material parameters
BEGIN_FECORE_CLASS(FEEFDMooneyRivlin, FEUncoupledMaterial)
	ADD_PARAMETER(m_MR.m_c1, "c1");
	ADD_PARAMETER(m_MR.m_c2, "c2");
	ADD_PARAMETER(m_EFD.m_beta, 3, "beta");
	ADD_PARAMETER(m_EFD.m_ksi , 3, "ksi" );
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEEFDMooneyRivlin::FEEFDMooneyRivlin(FEModel* pfem) : FEUncoupledMaterial(pfem), m_EFD(pfem), m_MR(pfem)
{
	// although we don't use K of the child materials, we need to set it to a non-zero value
	// otherwise FEBio will complain
	m_MR.m_K = 1.0;
	m_EFD.m_K = 1.0;
}

//-----------------------------------------------------------------------------
bool FEEFDMooneyRivlin::Init()
{
	if (FEUncoupledMaterial::Init() == false) return false;
	if (m_MR.Init() == false) return false;
	if (m_EFD.Init() == false) return false;
	return true;
}

//-----------------------------------------------------------------------------
void FEEFDMooneyRivlin::Serialize(DumpStream& ar)
{
	FEUncoupledMaterial::Serialize(ar);
	m_MR.Serialize(ar);
	m_EFD.Serialize(ar);
}

//-----------------------------------------------------------------------------
mat3ds FEEFDMooneyRivlin::DevStress(FEMaterialPoint& pt)
{
	return m_MR.DevStress(pt) + m_EFD.DevStress(pt);
}

//-----------------------------------------------------------------------------
tens4ds FEEFDMooneyRivlin::DevTangent(FEMaterialPoint& pt)
{
	return m_MR.DevTangent(pt) + m_EFD.DevTangent(pt);
}

//-----------------------------------------------------------------------------
//! calculate deviatoric strain energy density
double FEEFDMooneyRivlin::DevStrainEnergyDensity(FEMaterialPoint& pt)
{
    return m_MR.DevStrainEnergyDensity(pt) + m_EFD.DevStrainEnergyDensity(pt);
}
