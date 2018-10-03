#pragma once
#include <FECore/FEMaterial.h>

//-----------------------------------------------------------------------------
//! A material class describing the active fiber contraction
class FEActiveFiberContraction : public FEMaterial
{
public:
	FEActiveFiberContraction(FEModel* pfem);

	//! initialization
	bool Init() override;

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

	DECLARE_FECORE_CLASS();
};
