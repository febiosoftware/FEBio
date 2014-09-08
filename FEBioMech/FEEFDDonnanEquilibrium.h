#pragma once
#include "FECore/FEMaterial.h"
#include "FEDonnanEquilibrium.h"
#include "FEEllipsoidalFiberDistribution.h"

//-----------------------------------------------------------------------------
//! This class implements a material that consists of a continuous ellipsoidal fiber distribution
//! superposed on a charged (swelling) gel described by the equations of Donnan equilibrium

//! This material is orignally due to Gerard Ateshian and is used to model
//! articular cartilage.

class FEEFDDonnanEquilibrium : public FEElasticMaterial
{
public:
	FEEFDDonnanEquilibrium(FEModel* pfem) : FEElasticMaterial(pfem), m_Fib(pfem), m_DEQ(pfem) {}
	
public:
	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt);
		
	//! calculate tangent stiffness at material point
	virtual tens4ds Tangent(FEMaterialPoint& pt);
		
	//! calculate strain energy density at material point
	virtual double StrainEnergyDensity(FEMaterialPoint& pt);
    
	//! material parameter intialization and checking
	void Init();

	//! serialization
	void Serialize(DumpFile& ar);
		
public:
		
//	FEEllipsoidalFiberDistribution		m_Fib;
	FEEllipsoidalFiberDistributionOld	m_Fib;
	FEDonnanEquilibrium					m_DEQ;
		
	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
