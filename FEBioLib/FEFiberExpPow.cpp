#include "stdafx.h"
#include "FEFiberExpPow.h"

// define the material parameters
BEGIN_PARAMETER_LIST(FEFiberExpPow, FEElasticMaterial)
	ADD_PARAMETER(m_alpha, FE_PARAM_DOUBLE, "alpha");
	ADD_PARAMETER(m_beta, FE_PARAM_DOUBLE, "beta");
	ADD_PARAMETER(m_ksi , FE_PARAM_DOUBLE, "ksi" );
	ADD_PARAMETER(m_thd, FE_PARAM_DOUBLE, "theta");
	ADD_PARAMETER(m_phd, FE_PARAM_DOUBLE, "phi");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
// FEFiberExpPow
//-----------------------------------------------------------------------------

void FEFiberExpPow::Init()
{
	if (m_unstable) throw MaterialError("This fibrous material is unstable (collapses on itself) when used alone.  Combine it in a solid mixture with a material that can serve as a ground matrix.");
	if (m_ksi < 0) throw MaterialError("ksi must be positive.");
	if (m_beta < 2) throw MaterialError("beta must be >= 2.");
	if (m_alpha < 0) throw MaterialError("alpha must be >= 0.");

	// convert angles from degrees to radians
	double pi = 4*atan(1.0);
	double the = m_thd*pi/180.;
	double phi = m_phd*pi/180.;
	// fiber direction in local coordinate system (reference configuration)
	m_n0.x = cos(the)*sin(phi);
	m_n0.y = sin(the)*sin(phi);
	m_n0.z = cos(phi);
}

//-----------------------------------------------------------------------------
mat3ds FEFiberExpPow::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// deformation gradient
	mat3d &F = pt.F;
	double J = pt.J;
	
	// loop over all integration points
	vec3d n0, nt;
	double In_1, Wl;
	const double eps = 0;
	mat3ds C = pt.RightCauchyGreen();
	mat3ds s;
	
	// evaluate fiber direction in global coordinate system
	n0 = pt.Q*m_n0;
	
	// Calculate In = n0*C*n0
	In_1 = n0*(C*n0) - 1.0;
	
	// only take fibers in tension into consideration
	if (In_1 > eps)
	{
		// get the global spatial fiber direction in current configuration
		nt = F*n0;
		
		// calculate the outer product of nt
		mat3ds N = dyad(nt);
		
		// calculate strain energy derivative
		Wl = m_ksi*pow(In_1, m_beta-1.0)*exp(m_alpha*pow(In_1, m_beta));
		
		// calculate the fiber stress
		s = N*(2.0*Wl/J);
	}
	else
	{
		s.zero();
	}
	
	return s;
}

//-----------------------------------------------------------------------------
tens4ds FEFiberExpPow::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// deformation gradient
	mat3d &F = pt.F;
	double J = pt.J;
	
	// loop over all integration points
	vec3d n0, nt;
	double In_1, Wll;
	const double eps = 0;
	mat3ds C = pt.RightCauchyGreen();
	tens4ds c;
	
	// evaluate fiber direction in global coordinate system
	n0 = pt.Q*m_n0;
	
	// Calculate In = n0*C*n0
	In_1 = n0*(C*n0) - 1.0;
	
	// only take fibers in tension into consideration
	if (In_1 > eps)
	{
		// get the global spatial fiber direction in current configuration
		nt = F*n0;
		
		// calculate the outer product of nt
		mat3ds N = dyad(nt);
		tens4ds NxN = dyad1s(N);
		
		// calculate strain energy 2nd derivative
		double tmp = m_alpha*pow(In_1, m_beta);
		Wll = m_ksi*pow(In_1, m_beta-2.0)*((tmp+1)*m_beta-1.0)*exp(tmp);
		
		// calculate the fiber tangent
		c = NxN*(4.0*Wll/J);
	}
	else
	{
		c.zero();
	}
	
	return c;
}
