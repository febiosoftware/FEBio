#pragma once
#include "FEMaterial.h"

//-----------------------------------------------------------------------------
//! This class implements a material that consists of a continuous fiber distribution

//! This material is orignally due to Gerard Ateshian and is used to model
//! articular cartilage. The only difference is that it uses a Mooney-Rivlin matrix.

class FERandomFiberMooneyRivlin :	public FEIncompressibleMaterial
{
public:
	FERandomFiberMooneyRivlin(void);

public:
	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt);

	//! calculate tangent stiffness at material point
	virtual void Tangent(double D[6][6], FEMaterialPoint& pt);

	//! material parameter intialization and checking
	void Init();

public:
	double	m_c1;	// Mooney-Rivlin coefficient 1
	double	m_c2;	// Mooney-Rivlin coefficient 2

	double	m_beta[3];
	double	m_ksi[3];

	//--- active contraction stuff ---
	double	m_a[3];
	double	m_ac;
	// -------------------------------

	static	int	m_nres;	// integration rule

	static double	m_cth[];
	static double	m_sth[];
	static double	m_cph[];
	static double	m_sph[];

	// declare as registered
	DECLARE_REGISTERED(FERandomFiberMooneyRivlin);

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
