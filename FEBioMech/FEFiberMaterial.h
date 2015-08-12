#pragma once
#include "FECore/FEMaterial.h"

//-----------------------------------------------------------------------------
//! A material class describing the active fiber contraction
class FEActiveFiberContraction : public FEMaterial
{
public:
	FEActiveFiberContraction(FEModel* pfem);

	//! initialization
	void Init();

	//! calculate the fiber stress
	mat3ds FiberStress(FEMaterialPoint& mp);

	//! active contraction stiffness contribution
	tens4ds FiberStiffness(FEMaterialPoint& mp);

protected:
	double	m_ascl;		//!< activation scale factor
	double	m_Tmax;		//!< activation scale factor
	double	m_ca0;		//!< intracellular calcium concentration
	double	m_camax;	//!< peak calcium concentration
	double	m_beta;		//!< shape of peak isometric tension-sarcomere length relation
	double	m_l0;		//!< unloaded length
	double	m_refl;		//!< sarcomere length

	DECLARE_PARAMETER_LIST();
};

//-----------------------------------------------------------------------------
//! Base class for fiber materials.
class FEFiberMaterial
{
public:
	//! Constructor
	FEFiberMaterial();

	//! Calculate the fiber stress
	mat3ds Stress(FEMaterialPoint& mp);

	//! Calculate the fiber tangent
	tens4ds Tangent(FEMaterialPoint& mp);

	//! Calculate the fiber strain energy density
	double StrainEnergyDensity(FEMaterialPoint& mp);
    
public:
	double	m_c3;		//!< Exponential stress coefficient
	double	m_c4;		//!< Fiber uncrimping coefficient
	double	m_c5;		//!< Modulus of straightened fibers

	double	m_lam1;		//!< fiber stretch for straightened fibers
};
