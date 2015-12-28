#pragma once
#include "FEElasticMaterial.h"

// The following file contains the integration points and weights
// for the integration over a unit sphere in spherical coordinates
#include "geodesic.h"

class FEFiberNeoHookean : public FEElasticMaterial
{
public:
	FEFiberNeoHookean(FEModel* pfem);

public:
	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt);

	//! calculate tangent stiffness at material point
	virtual tens4ds Tangent(FEMaterialPoint& pt);

	//! calculate strain energy density at material point
	virtual double StrainEnergyDensity(FEMaterialPoint& mp);
    
	//! material parameter intialization and checking
	bool Init();

public:
	double	m_E;	//<! Young's modulus
	double	m_v;	//<! Poisson's ratio

	//--- active contraction stuff ---
	double	m_a[3];
	double	m_ac;
	// -------------------------------

	// numerical quadrature stuff
	int     m_nres;	// integration rule

	double	m_cth[NSTH];
	double	m_sth[NSTH];
	double	m_cph[NSTH];
	double	m_sph[NSTH];
    double  m_w[NSTH];
    bool    m_bfirst;

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
