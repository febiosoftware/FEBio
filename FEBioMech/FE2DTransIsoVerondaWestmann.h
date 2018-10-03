#pragma once
#include "FEUncoupledMaterial.h"

//-----------------------------------------------------------------------------
//! 2D transversely isotropic Veronda-Westmann

//! This class describes a transversely isotropic matrix where the base material
//! is Veronda-Westmann. The difference between this material and the FETransIsoVerondaWestmann
//! material is that in this material the fibers lie in the plane that is perpendicular
//! to the transverse axis. 

class FE2DTransIsoVerondaWestmann :	public FEUncoupledMaterial
{
	enum { NSTEPS = 12 };	// nr of integration steps

public:
	// material parameters
	double	m_c1;	//!< Veronda-Westmann parameter c1
	double	m_c2;	//!< Veronda-Westmann parameter c2

	// fiber parameters
	double	m_c3;
	double	m_c4;
	double	m_c5;
	double	m_lam1;

public:
	//! constructor
	FE2DTransIsoVerondaWestmann(FEModel* pfem);
	
	//! calculate deviatoric stress at material point
	virtual mat3ds DevStress(FEMaterialPoint& pt) override;

	//! calculate deviatoric tangent stiffness at material point
	virtual tens4ds DevTangent(FEMaterialPoint& pt) override;

	// declare parameter list
	DECLARE_FECORE_CLASS();

protected:
	static double	m_cth[NSTEPS];
	static double	m_sth[NSTEPS];

	double	m_w[2];
};
