#include "FEUncoupledMaterial.h"

//-----------------------------------------------------------------------------
//! Material class for single fiber, tension only
//! Exponential-power law

class FEFiberExpPowUncoupled : public FEUncoupledMaterial
{
public:
	FEFiberExpPowUncoupled(FEModel* pfem) : FEUncoupledMaterial(pfem) { m_thd = 0; m_phd = 90; }
	
	//! Initialization
	void Init();

	//! Cauchy stress
	virtual mat3ds DevStress(FEMaterialPoint& mp);
	
	// Spatial tangent
	virtual tens4ds DevTangent(FEMaterialPoint& mp);
	
	// declare the parameter list
	DECLARE_PARAMETER_LIST();
	
public:
	double	m_alpha;	// coefficient of (In-1) in exponential
	double	m_beta;		// power of (In-1) in exponential
	double	m_ksi;		// fiber modulus
	double	m_thd;		// theta angle for fiber orientation (local coordinates system)
	double	m_phd;		// phi angle for fiber orientation (local coordinates system)
	vec3d	m_n0;		// unit vector along fiber direction (local coordinate system)
};
