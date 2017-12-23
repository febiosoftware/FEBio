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
	bool Init() override;

	//! serialization
	void Serialize(DumpStream& ar) override;

public:
	//! Calculate the deviatoric stress
	mat3ds DevStress(FEMaterialPoint& pt) override;

	//! Calculate deviatoric tangent
	tens4ds DevTangent(FEMaterialPoint& pt) override;

	//! calculate deviatoric strain energy density
	double DevStrainEnergyDensity(FEMaterialPoint& mp) override;
    
public:
	FEMooneyRivlin	m_MR;
	FEEFDUncoupled	m_EFD;

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
