#pragma once
#include "FEMaterial.h"

class FE2DFiberNeoHookean :	public FEElasticMaterial
{
	enum { NSTEPS = 12 };	// nr of integration steps

public:
	FE2DFiberNeoHookean();

public:
	double	m_E;	//!< Young's modulus
	double	m_v;	//!< Poisson's ratio

	//--- active contraction stuff ---
	double	m_a[2];
	double	m_ac;
	// -------------------------------

public:
	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt);

	//! calculate tangent stiffness at material point
	virtual tens4ds Tangent(FEMaterialPoint& pt);

	//! data initialization and checking
	void Init();

	//! return bulk modulus
	double BulkModulus() { return m_E/(3.0*(1.0 - 2.0*m_v));}

	// declare as registered
	DECLARE_REGISTERED(FE2DFiberNeoHookean);

	// declare the parameter list
	DECLARE_PARAMETER_LIST();


protected:
	static double	m_cth[NSTEPS];
	static double	m_sth[NSTEPS];
};
