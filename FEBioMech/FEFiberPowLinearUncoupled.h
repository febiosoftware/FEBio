#pragma once
#include "FEElasticFiberMaterialUC.h"

//-----------------------------------------------------------------------------
//! Material class for single fiber, tension only
//! Power-law linear response (uncoupled)

class FEFiberPowLinearUncoupled : public FEElasticFiberMaterialUC
{
public:
    FEFiberPowLinearUncoupled(FEModel* pfem);

	// validation
	bool Validate() override;
    
    //! Cauchy stress
    virtual mat3ds DevStress(FEMaterialPoint& mp, const vec3d& n0) override;
    
    // Spatial tangent
    virtual tens4ds DevTangent(FEMaterialPoint& mp, const vec3d& n0) override;
    
    //! Strain energy density
    virtual double DevStrainEnergyDensity(FEMaterialPoint& mp, const vec3d& n0) override;
    
public:
    double	m_E;		// fiber modulus
    double  m_lam0;     // stretch ratio at end of toe region
	double  m_beta;     // power law exponent in toe region

private:
    double  m_ksi;      // power law coefficient in toe region
	double  m_I0;       // m_lam0^2
    double  m_b;        // coefficient in linear region

	// declare the parameter list
	DECLARE_FECORE_CLASS();
};
