#pragma once
#include "FECore/FEMaterial.h"
#include "FENeoHookean.h"
#include "FEEllipsoidalFiberDistribution.h"

class FEEFDNeoHookean :	public FEElasticMaterial
{
public:
	FEEFDNeoHookean(FEModel* pfem);

public:
	double	m_E;	//!< Young's modulus
	double	m_v;	//!< Poisson's ratio
	double	m_ksi[3];	//!< ksi
	double	m_beta[3];	//!< beta

public:
	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt) override;

	//! calculate tangent stiffness at material point
	virtual tens4ds Tangent(FEMaterialPoint& pt) override;

	//! calculate strain energy density at material point
	virtual double StrainEnergyDensity(FEMaterialPoint& pt) override;
    
	//! data initialization
	bool Init() override;

	//! serialization
	void Serialize(DumpStream& ar) override;

public:

	FEEllipsoidalFiberDistribution	m_EFD;
	FENeoHookean					m_NH;

	// declare the parameter list
	DECLARE_FECORE_CLASS();
};


//----------------------------------------------------------------
// "old" version using full-sphere integration
class FEEFDNeoHookeanOld :	public FEElasticMaterial
{
public:
	FEEFDNeoHookeanOld(FEModel* pfem) : FEElasticMaterial(pfem), m_EFD(pfem), m_NH(pfem) {}

public:
	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt) override;

	//! calculate tangent stiffness at material point
	virtual tens4ds Tangent(FEMaterialPoint& pt) override;

	//! calculate strain energy density at material point
	virtual double StrainEnergyDensity(FEMaterialPoint& pt) override;
    
	//! data initialization and checking
	bool Init() override;

	//! serialization
	void Serialize(DumpStream& ar) override;

public:

	FEEllipsoidalFiberDistributionOld	m_EFD;
	FENeoHookean						m_NH;

	// declare the parameter list
	DECLARE_FECORE_CLASS();
};
