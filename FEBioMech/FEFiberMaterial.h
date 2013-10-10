#pragma once
#include "FECore/FEMaterial.h"

//-----------------------------------------------------------------------------
//! Base class for fiber materials.
class FEFiberMaterial : public FEMaterial
{
public:
	//! Constructor
	FEFiberMaterial();

	//! Initialization
	void Init();

	//! Calculate the fiber stress
	mat3ds Stress(FEMaterialPoint& mp);

	//! Calculate the fiber tangent
	tens4ds Tangent(FEMaterialPoint& mp);

public:
	double	m_c3;		//!< Exponential stress coefficient
	double	m_c4;		//!< Fiber uncrimping coefficient
	double	m_c5;		//!< Modulus of straightened fibers

	double	m_lam1;		//!< fiber stretch for straightened fibers

	//--- time varying elastance active contraction data ---
	double	m_ascl;		//!< activation scale factor
	double	m_Tmax;		//!< activation scale factor
	double	m_ca0;		//!< intracellular calcium concentration
	double	m_camax;	//!< peak calcium concentration
	double	m_beta;		//!< shape of peak isometric tension-sarcomere length relation
	double	m_l0;		//!< unloaded length
	double	m_refl;		//!< sarcomere length

	//--- pre-strain data ---
	double	m_lcur;		//!< in-situ current stretch

	DECLARE_PARAMETER_LIST();
};
