#pragma once
#include "FEElasticMaterial.h"

//-----------------------------------------------------------------------------
class FE2DFiberNeoHookean :	public FEElasticMaterial
{
	enum { NSTEPS = 12 };	// nr of integration steps

public:
	FE2DFiberNeoHookean(FEModel* pfem);

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

	// declare the parameter list
	DECLARE_PARAMETER_LIST();

protected:
	static double	m_cth[NSTEPS];
	static double	m_sth[NSTEPS];
};
