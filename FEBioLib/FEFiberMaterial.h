#pragma once
#include "FECore/FEMaterial.h"

//-----------------------------------------------------------------------------
//! Base class for fiber materials.
class FEFiberMaterial : public FEMaterial
{
public:
	//! Constructor
	FEFiberMaterial()
	{
		m_ascl = 0;
		m_c3 = m_c4 = m_c5 = 0;
		m_lam1 = 1;
	}

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
	double	m_ca0;		//!< intracellular calcium concentration
	double	m_beta;		//!< shape of peak isometric tension-sarcomere length relation
	double	m_l0;		//!< unloaded length
	double	m_refl;		//!< sarcomere length

	DECLARE_PARAMETER_LIST();
};
