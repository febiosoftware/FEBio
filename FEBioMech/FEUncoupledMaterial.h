#pragma once
#include "FEElasticMaterial.h"

//-----------------------------------------------------------------------------
//! Base class for uncoupled hyperelastic material formulations.

//! In FEBio, for uncoupled materials it is assumed that the strain energy function 
//! is a sum of two terms, a deviatoric strain energy term, which only depends on the 
//! deviatoric right Cauchy-Green tensor C_tilde, and a volumetric term which only 
//! depends on J, the determinant of the deformation gradient. The total Cauchy stress 
//! is therefore also a sum of two contributions, namely the deviatoric Cauchy stress
//! and a pressure term where p = dU/dJ. 

//! One of the main motivations for the alternative material interface is that some
//! finite element implementations can take advantage of the uncoupling. For example,
//! the three-field formulation integrates the two different material terms differently
//! in order to avoid locking problems for (nearly-) incompressible materials. 

//! When implementing a new material derived from this base class, the developer needs
//! to provide the deviatoric stress function and the pressure function as well as their
//! derivatives. 

class FEUncoupledMaterial :	public FEElasticMaterial
{
public:
	//! constructor
	FEUncoupledMaterial(FEModel* pfem);

public:

//----------------->
	// The following four functions need to be overloaded
	// for each material derived from this class.

	//! Deviatoric Cauchy stress
	virtual mat3ds DevStress(FEMaterialPoint& mp) = 0;

	//! Deviatoric spatial Tangent
	virtual tens4ds DevTangent(FEMaterialPoint& mp) = 0;

	//! Deviatoric strain energy density
	virtual double DevStrainEnergyDensity(FEMaterialPoint& mp) { return 0; }
    
	//! strain energy density U(J)
	virtual double U(double J) { double lnJ = log(J); return 0.5*m_K*lnJ*lnJ; }
    
	//! pressure, i.e. first derivative of U(J)
	virtual double UJ(double J) { return m_K*log(J)/J; }

	//! second derivative of U(J) 
	virtual double UJJ(double J) { return m_K*(1-log(J))/(J*J); }

	// use these for NIKE3D's Ogden material
//	double U  (double J) { return 0.25*m_K*(J*J - 2.0*log(J) - 1.0); }
//	double UJ (double J) { return 0.5*m_K*(J - 1.0/J); }
//	double UJJ(double J) { return 0.5*m_K*(1 + 1.0/(J*J)); }

	// Use these to obtain similar results than ABAQUS
//	double U  (double J) { return 0.5*m_K*(J-1)*(J-1); }
//	double UJ (double J) { return m_K*(J-1); }
//	double UJJ(double J) { return m_K; }
//----------------->

	// incompressibility constraint fnc and derivs
	double h  (double J) { return log(J); }
	double hp (double J) { return 1.0 / J; }
	double hpp(double J) { return -1.0 / (J*J); }

public:
	//! total Cauchy stress (do not overload!)
	mat3ds Stress(FEMaterialPoint& mp);

	//! total spatial tangent (do not overload!)
	tens4ds Tangent(FEMaterialPoint& mp);

	//! calculate strain energy (do not overload!)
	double StrainEnergyDensity(FEMaterialPoint& pt);
    
	//! material initialization
	void Init();

public:
	bool	m_blaugon;	//!< augmented lagrangian flag
	double	m_atol;		//!< augmented lagrangian tolerance
	double	m_K;		//!< bulk modulus

	DECLARE_PARAMETER_LIST();
};
