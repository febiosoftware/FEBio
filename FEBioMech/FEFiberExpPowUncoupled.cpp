#include "stdafx.h"
#include "FEFiberExpPowUncoupled.h"

// define the material parameters
BEGIN_PARAMETER_LIST(FEFiberExpPowUncoupled, FEUncoupledMaterial)
	ADD_PARAMETER(m_alpha, FE_PARAM_DOUBLE, "alpha");
	ADD_PARAMETER(m_beta, FE_PARAM_DOUBLE, "beta");
	ADD_PARAMETER(m_ksi , FE_PARAM_DOUBLE, "ksi" );
	ADD_PARAMETER(m_thd, FE_PARAM_DOUBLE, "theta");
	ADD_PARAMETER(m_phd, FE_PARAM_DOUBLE, "phi");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
// FEFiberExpPowUncoupled
//-----------------------------------------------------------------------------

void FEFiberExpPowUncoupled::Init()
{
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
mat3ds FEFiberExpPowUncoupled::DevStress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// deformation gradient
	double J = pt.m_J;
	mat3d F = pt.m_F*pow(J,-1.0/3.0);
	
	// loop over all integration points
	vec3d n0, nt;
	double In_1, Wl;
	const double eps = 0;
	mat3ds C = pt.DevRightCauchyGreen();
	mat3ds s;
	
	// evaluate fiber direction in global coordinate system
	n0 = pt.m_Q*m_n0;
	
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
	
	return s.dev();
}

//-----------------------------------------------------------------------------
tens4ds FEFiberExpPowUncoupled::DevTangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// deformation gradient
	double J = pt.m_J;
	mat3d F = pt.m_F*pow(J,-1.0/3.0);
	
	// loop over all integration points
	vec3d n0, nt;
	double In_1, Wl, Wll;
	const double eps = 0;
	mat3ds C = pt.DevRightCauchyGreen();
	mat3ds s;
	tens4ds c;
	
	// evaluate fiber direction in global coordinate system
	n0 = pt.m_Q*m_n0;
	
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
		
		// calculate strain energy derivatives
		double tmp = m_alpha*pow(In_1, m_beta);
		Wl = m_ksi*pow(In_1, m_beta-1.0)*exp(m_alpha*pow(In_1, m_beta));
		Wll = m_ksi*pow(In_1, m_beta-2.0)*((tmp+1)*m_beta-1.0)*exp(tmp);
		
		// calculate the fiber stress
		s = N*(2.0*Wl/J);
		
		// calculate the fiber tangent
		c = NxN*(4.0*Wll/J);
		
		// This is the final value of the elasticity tensor
		mat3dd I(1);
		tens4ds IxI = dyad1s(I);
		tens4ds I4  = dyad4s(I);
		c += ((I4+IxI/3.0)*s.tr() - dyad1s(I,s))*(2./3.)
		- (ddots(IxI, c)-IxI*(c.tr()/3.))/3.;
	}
	else
	{
		c.zero();
	}
	
	return c;
}

//-----------------------------------------------------------------------------
double FEFiberExpPowUncoupled::DevStrainEnergyDensity(FEMaterialPoint& mp)
{
    double sed = 0.0;
    
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// loop over all integration points
	vec3d n0;
	double In_1;
	const double eps = 0;
	mat3ds C = pt.DevRightCauchyGreen();
	
	// evaluate fiber direction in global coordinate system
	n0 = pt.m_Q*m_n0;
	
	// Calculate In = n0*C*n0
	In_1 = n0*(C*n0) - 1.0;
	
	// only take fibers in tension into consideration
	if (In_1 > eps)
	{
		// calculate strain energy derivative
        if (m_alpha > 0) {
            sed = m_ksi/(m_alpha*m_beta)*(exp(m_alpha*pow(In_1, m_beta))-1);
        }
        else
            sed = m_ksi/m_beta*pow(In_1, m_beta);
	}
    
    return sed;
}
