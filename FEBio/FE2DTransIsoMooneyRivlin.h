#pragma once
#include "FEMaterial.h"

//-----------------------------------------------------------------------------
//! 2D transversely isotropic Mooney-Rivlin

//! This class describes a transversely isotropic matrix where the base material
//! is Mooney-Rivlin. The difference between this material and the FETransIsoMooneyRivlin
//! material is that in this material the fibers lie in the plane that is perpendicular
//! to the transverse axis. 

class FE2DTransIsoMooneyRivlin : public FETransverselyIsotropic
{
	enum { NSTEPS = 12 };	// nr of integration steps

public:
	// material parameters
	double	m_c1;	//!< Mooney-Rivlin parameter c1
	double	m_c2;	//!< Mooney-Rivlin parameter c2

public:
	//! constructor
	FE2DTransIsoMooneyRivlin();
	
	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt);

	//! calculate tangent stiffness at material point
	virtual void Tangent(double D[6][6], FEMaterialPoint& pt);

	// declare as registered
	DECLARE_REGISTERED(FE2DTransIsoMooneyRivlin);

	// declare parameter list
	DECLARE_PARAMETER_LIST();

protected:
	static double	m_cth[NSTEPS];
	static double	m_sth[NSTEPS];
};
