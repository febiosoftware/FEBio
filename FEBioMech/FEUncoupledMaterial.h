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
    virtual double U(double J) {
        switch (m_npmodel) {
            case 0: return 0.5*m_K*pow(log(J),2); break;    // FEBio default
            case 1: return 0.25*m_K*(J*J - 2.0*log(J) - 1.0); break;    // NIKE3D's Ogden material
            case 2: return 0.5*m_K*(J-1)*(J-1); break;      // ABAQUS
            default: { assert(false); return 0; }
        }
    }
	//! pressure, i.e. first derivative of U(J)
	virtual double UJ(double J) {
        switch (m_npmodel) {
            case 0: return m_K*log(J)/J; break;
            case 1: return 0.5*m_K*(J - 1.0/J); break;
            case 2: return m_K*(J-1); break;
			default: { assert(false); return 0; }
		}
    }

	//! second derivative of U(J) 
	virtual double UJJ(double J) {
        switch (m_npmodel) {
            case 0: return m_K*(1-log(J))/(J*J); break;
            case 1: return 0.5*m_K*(1 + 1.0/(J*J)); break;
            case 2: return m_K; break;
			default: { assert(false); return 0; }
		}
    }

	// incompressibility constraint fnc and derivs
	double h  (double J) { return log(J); }
	double hp (double J) { return 1.0 / J; }
	double hpp(double J) { return -1.0 / (J*J); }

public:
	//! total Cauchy stress (do not overload!)
	mat3ds Stress(FEMaterialPoint& mp) override;

	//! total spatial tangent (do not overload!)
	tens4ds Tangent(FEMaterialPoint& mp) override;

	//! calculate strain energy (do not overload!)
	double StrainEnergyDensity(FEMaterialPoint& pt) override;

	// Create material point data
	FEMaterialPoint* CreateMaterialPointData() override;
    
public:
	double	m_K;		//!< bulk modulus

	bool	m_blaugon;		//!< augmented lagrangian flag
	double	m_augtol;		//!< augmented lagrangian tolerance
	int		m_naugmin;		//!< minimum number of augmentations
	int		m_naugmax;		//!< max number of augmentations
    int     m_npmodel;      //!< pressure model for U(J)

	DECLARE_FECORE_CLASS();
};
