#pragma once
#include "FEUncoupledMaterial.h"
#include "FEVerondaWestmann.h"
#include "FEEFDUncoupled.h"

//-----------------------------------------------------------------------------
//! This class implements a material that consists of a Veronda-Westmann matrix and
//! a continuous EFD fiber distribution.
class FEEFDVerondaWestmann : public FEUncoupledMaterial
{
public:
	// constructor
	FEEFDVerondaWestmann(FEModel* pfem);

	//! material initialization
	bool Init();

	//! deviatoric stress
	mat3ds DevStress(FEMaterialPoint& pt);

	//! deviatoric tangent
	tens4ds DevTangent(FEMaterialPoint& pt);

	//! calculate deviatoric strain energy density
	double DevStrainEnergyDensity(FEMaterialPoint& mp);
    
public:
	FEVerondaWestmann	m_VW;
	FEEFDUncoupled		m_EFD;

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
