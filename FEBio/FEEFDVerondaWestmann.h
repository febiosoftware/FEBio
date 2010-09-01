#pragma once
#include "FEUncoupledMaterial.h"
#include "FEVerondaWestmann.h"
#include "FEEllipsoidalFiberDistribution.h"

//-----------------------------------------------------------------------------
//! This class implements a material that consists of a continuous fiber distribution

//! This material is orignally due to Gerard Ateshian and is used to model
//! articular cartilage. The only difference is that it uses a Veronda-Westmann matrix.

class FEEFDVerondaWestmann:	public FEUncoupledMaterial
{
public:
	FEEFDVerondaWestmann() {}

public:
	//! calculate deviatoric stress at material point
	virtual mat3ds DevStress(FEMaterialPoint& pt);

	//! calculate deviatoric tangent stiffness at material point
	virtual tens4ds DevTangent(FEMaterialPoint& pt);

	//! material parameter intialization and checking
	void Init();

	//! return bulk modulus
	double BulkModulus();
	
public:
	FEEllipsoidalFiberDistribution	m_EFD;
	FEVerondaWestmann				m_VW;

	// declare as registered
	DECLARE_REGISTERED(FEEFDVerondaWestmann);

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
