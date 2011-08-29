#pragma once
#include "FECore/FEElasticMaterial.h"

class FEFiberNeoHookean : public FEElasticMaterial
{
public:
	FEFiberNeoHookean();

public:
	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt);

	//! calculate tangent stiffness at material point
	virtual tens4ds Tangent(FEMaterialPoint& pt);

	//! material parameter intialization and checking
	void Init();

	//! return bulk modulus
	double BulkModulus() { return m_E/(3.0*(1.0 - 2.0*m_v));}

public:
	double	m_E;	//<! Young's modulus
	double	m_v;	//<! Poisson's ratio

	//--- active contraction stuff ---
	double	m_a[3];
	double	m_ac;
	// -------------------------------

	// numerical quadrature stuff
	static	int	m_nres;	// integration rule

	static double	m_cth[];
	static double	m_sth[];
	static double	m_cph[];
	static double	m_sph[];

	// declare as registered
	DECLARE_REGISTERED(FEFiberNeoHookean);

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
