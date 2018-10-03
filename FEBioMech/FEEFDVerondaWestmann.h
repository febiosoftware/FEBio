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
	bool Init() override;

	//! serialization
	void Serialize(DumpStream& ar) override;

	//! deviatoric stress
	mat3ds DevStress(FEMaterialPoint& pt) override;

	//! deviatoric tangent
	tens4ds DevTangent(FEMaterialPoint& pt) override;

	//! calculate deviatoric strain energy density
	double DevStrainEnergyDensity(FEMaterialPoint& mp) override;
    
public:
	FEVerondaWestmann	m_VW;
	FEEFDUncoupled		m_EFD;

	// declare the parameter list
	DECLARE_FECORE_CLASS();
};
