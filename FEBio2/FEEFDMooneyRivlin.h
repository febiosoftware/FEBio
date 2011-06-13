#pragma once
#include "FEBioLib/FEUncoupledMaterial.h"

//-----------------------------------------------------------------------------
//! This class implements a material that consists of a continuous fiber distribution

//! This material is orignally due to Gerard Ateshian and is used to model
//! articular cartilage. The only difference is that it uses a Mooney-Rivlin matrix.

class FEEFDMooneyRivlin :	public FEUncoupledMaterial
{
public:
	FEEFDMooneyRivlin(void);

public:
	//! calculate deviatoric stress at material point
	mat3ds DevStress(FEMaterialPoint& pt);

	//! calculate deviatoric tangent stiffness at material point
	tens4ds DevTangent(FEMaterialPoint& pt);

	//! material parameter intialization and checking
	void Init();

protected:
	//! Calculate (deviatoric) fiber stress
	mat3ds FiberStress(FEMaterialPoint& mp);

	//! Calculate (deviatoric) fiber tangent
	tens4ds FiberTangent(FEMaterialPoint& mp);

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
	static double	m_w[];

	// declare as registered
	DECLARE_REGISTERED(FEEFDMooneyRivlin);

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
