#pragma once
#include "FEUncoupledMaterial.h"
#include "FEMooneyRivlin.h"
#include "FEEFDUncoupled.h"

//-----------------------------------------------------------------------------
//! This class implements a material that consists of a Mooney-Rivlin matrix and
//! a continuous EFD fiber distribution.
class FEEFDMooneyRivlin : public FEUncoupledMaterial
{
public:
	// constructor
	FEEFDMooneyRivlin(FEModel* pfem);
	
	//! Data initialization
	bool Init();

	//! serialization
	void Serialize(DumpStream& ar);

public:
	//! Calculate the deviatoric stress
	mat3ds DevStress(FEMaterialPoint& pt);

	//! Calculate deviatoric tangent
	tens4ds DevTangent(FEMaterialPoint& pt);

	//! calculate deviatoric strain energy density
	double DevStrainEnergyDensity(FEMaterialPoint& mp);
    
public:
	FEMooneyRivlin	m_MR;
	FEEFDUncoupled	m_EFD;

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
