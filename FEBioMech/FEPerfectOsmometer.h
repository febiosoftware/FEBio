#pragma once
#include "FEElasticMaterial.h"

//-----------------------------------------------------------------------------
//! Material class that implements the equilibrium of a perfect osmometer.
//! When used on its own (not in a solid mixture), this materials
//! is intrinsically unstable
class FEPerfectOsmometer : public FEElasticMaterial
{
public:
	// constructor
	FEPerfectOsmometer(FEModel* pfem);
	
	//! Initialization routine
	bool Init() override;

	//! Returns the Cauchy stress
	virtual mat3ds Stress(FEMaterialPoint& mp) override;

	//! Returs the spatial tangent
	virtual tens4ds Tangent(FEMaterialPoint& mp) override;

	// declare the parameter list
	DECLARE_FECORE_CLASS();
	
public:
	double	m_phiwr;	//!< fluid volume fraction in reference configuration
	double	m_iosm;		//!< internal osmolarity in reference configuration
	double	m_Rgas;		//!< universal gas constant
	double	m_Tabs;		//!< absolute temperature
	double	m_bosm;		//!< bath osmolarity
};
