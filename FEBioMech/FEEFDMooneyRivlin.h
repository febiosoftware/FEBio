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
	void Init();

public:
	//! Calculate the deviatoric stress
	mat3ds DevStress(FEMaterialPoint& pt);

	//! Calculate deviatoric tangent
	tens4ds DevTangent(FEMaterialPoint& pt);

public:
	FEMooneyRivlin	m_MR;
	FEEFDUncoupled	m_EFD;

	// declare as registered
	DECLARE_REGISTERED(FEEFDMooneyRivlin);

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
