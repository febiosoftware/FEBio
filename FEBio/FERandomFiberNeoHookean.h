#pragma once
#include "FEMaterial.h"

class FERandomFiberNeoHookean :	public FEElasticMaterial
{
public:
	FERandomFiberNeoHookean() {}

public:
	double	m_E;	//!< Young's modulus
	double	m_v;	//!< Poisson's ratio
	double	m_ksi[3];	//!< ksi
	double	m_beta[3];	//!< beta

public:
	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt);

	//! calculate tangent stiffness at material point
	virtual void Tangent(double D[6][6], FEMaterialPoint& pt);

	//! data initialization and checking
	void Init();

	//! return bulk modulus
	double BulkModulus() { return m_E/(3.0*(1.0 - 2.0*m_v));}

public:
	static	int	m_nres;	// integration rule
	static double	m_cth[];
	static double	m_sth[];
	static double	m_cph[];
	static double	m_sph[];


	// declare as registered
	DECLARE_REGISTERED(FERandomFiberNeoHookean);

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
