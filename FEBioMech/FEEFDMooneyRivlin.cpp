#include "stdafx.h"
#include "FEEFDMooneyRivlin.h"

// define the material parameters
BEGIN_PARAMETER_LIST(FEEFDMooneyRivlin, FEUncoupledMaterial)
	ADD_PARAMETER(m_MR.c1, FE_PARAM_DOUBLE, "c1");
	ADD_PARAMETER(m_MR.c2, FE_PARAM_DOUBLE, "c2");
	ADD_PARAMETERV(m_EFD.m_beta, FE_PARAM_DOUBLE, 3, "beta");
	ADD_PARAMETERV(m_EFD.m_ksi , FE_PARAM_DOUBLE, 3, "ksi" );
END_PARAMETER_LIST();

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
