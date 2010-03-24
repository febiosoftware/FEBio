#pragma once
#include "FEMaterial.h"
#include "FEDonnanEquilibrium.h"
#include "FEEllipsoidalFiberDistribution.h"

//-----------------------------------------------------------------------------
//! This class implements a material that consists of a continuous ellipsoidal fiber distribution
//! superposed on a charged (swelling) gel described by the equations of Donnan equilibrium

//! This material is orignally due to Gerard Ateshian and is used to model
//! articular cartilage.

class FERandomFiberDonnanEquilibrium :	public FEElasticMaterial
{
public:
	FERandomFiberDonnanEquilibrium(void);
	
public:
	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt);
		
	//! calculate tangent stiffness at material point
	virtual tens4ds Tangent(FEMaterialPoint& pt);
		
	//! material parameter intialization and checking
	void Init();

	//! return bulk modulus
	virtual double BulkModulus();
		
public:
		
	FEEllipsoidalFiberDistribution	m_Fib;
	FEDonnanEquilibrium				m_DEQ;
		
	// declare as registered
	DECLARE_REGISTERED(FERandomFiberDonnanEquilibrium);
		
	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
