#include "FEElasticFiberMaterialUC.h"
#include "FEFiberMaterialPoint.h"

BEGIN_PARAMETER_LIST(FEElasticFiberMaterialUC, FEUncoupledMaterial)
	ADD_PARAMETER(m_thd, FE_PARAM_DOUBLE, "theta");
	ADD_PARAMETER(m_phd, FE_PARAM_DOUBLE, "phi");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEElasticFiberMaterialUC::FEElasticFiberMaterialUC(FEModel* pfem) : FEUncoupledMaterial(pfem) 
{
	m_thd = 0.0;
	m_phd = 90.0;
}

//-----------------------------------------------------------------------------
FEMaterialPoint* FEElasticFiberMaterialUC::CreateMaterialPointData()
{
	FEFiberMaterialPoint* fp = new FEFiberMaterialPoint(FEUncoupledMaterial::CreateMaterialPointData());

	// Some fiber materials defined the theta,phi parameters for setting the fiber vector
	// Although this is deprecated, we still support it here for backward compatibility
	if ((m_thd != 0.0) || (m_phd != 90.0))
	{
		// convert angles from degrees to radians
		double pi = 4 * atan(1.0);
		double the = m_thd*pi / 180.;
		double phi = m_phd*pi / 180.;

		// fiber direction in local coordinate system (reference configuration)
		vec3d n0;
		n0.x = cos(the)*sin(phi);
		n0.y = sin(the)*sin(phi);
		n0.z = cos(phi);
		n0.unit();
		fp->m_n0 = n0;
	}

	return fp;
}

//-----------------------------------------------------------------------------
vec3d FEElasticFiberMaterialUC::GetFiberVector(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	FEFiberMaterialPoint& fp = *mp.ExtractData<FEFiberMaterialPoint>();

	return pt.m_Q*fp.m_n0;
}

//-----------------------------------------------------------------------------
// FEFiberExponentialPowerUC
//-----------------------------------------------------------------------------

// define the material parameters
BEGIN_PARAMETER_LIST(FEFiberExponentialPowerUC, FEElasticFiberMaterialUC)
    ADD_PARAMETER2(m_alpha, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "alpha");
    ADD_PARAMETER2(m_beta, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(2.0), "beta");
    ADD_PARAMETER(m_ksi , FE_PARAM_DOUBLE, "ksi" );
    ADD_PARAMETER2(m_mu   , FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "mu"   );
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEFiberExponentialPowerUC::FEFiberExponentialPowerUC(FEModel* pfem) : FEElasticFiberMaterialUC(pfem)
{
	m_alpha = 0; 
	m_beta = 2; 
	m_ksi = 0; 
	m_mu = 0;
}

//-----------------------------------------------------------------------------
bool FEFiberExponentialPowerUC::Validate()
{
	if ((4*m_ksi + 2*m_mu) < 0) return MaterialError("4*ksi+2*mu must be positive.");
    return FEElasticFiberMaterialUC::Validate();
}

//-----------------------------------------------------------------------------
mat3ds FEFiberExponentialPowerUC::DevStress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	double J = pt.m_J;
	// distortional part of deformation gradient
	mat3d F = pt.m_F*pow(J,-1.0/3.0);
	
	// loop over all integration points
	mat3ds C = pt.DevRightCauchyGreen();
	mat3ds s;
	
	// fiber direction in global coordinate system
	vec3d n0 = GetFiberVector(mp);
	
	// Calculate In = n0*C*n0
	double In_1 = n0*(C*n0) - 1.0;
	
	// only take fibers in tension into consideration
	const double eps = 0;
	if (In_1 > eps)
	{
		// get the global spatial fiber direction in current configuration
		vec3d nt = F*n0;
		
		// calculate the outer product of nt
		mat3ds N = dyad(nt);
		
		// calculate strain energy derivative
		double Wl = m_ksi*pow(In_1, m_beta-1.0)*exp(m_alpha*pow(In_1, m_beta));
		
		// calculate the fiber stress
		s = N*(2.0*Wl/J);
        
        // add the contribution from shear
        mat3ds BmI = pt.DevLeftCauchyGreen() - mat3dd(1);
        s += (N*BmI + BmI*N)*(m_mu/2/J);
	}
	else
	{
		s.zero();
	}
	
	return s.dev();
}

//-----------------------------------------------------------------------------
tens4ds FEFiberExponentialPowerUC::DevTangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	double J = pt.m_J;
	// distortional part of deformation gradient
	mat3d F = pt.m_F*pow(J,-1.0/3.0);
	
	const double eps = 0;
	mat3ds C = pt.DevRightCauchyGreen();
    mat3ds s;
	tens4ds c;
	
	// fiber direction in global coordinate system
	vec3d n0 = GetFiberVector(mp);
	
	// Calculate In = n0*C*n0
	double In_1 = n0*(C*n0) - 1.0;
	
	// only take fibers in tension into consideration
	if (In_1 > eps)
	{
		// get the global spatial fiber direction in current configuration
		vec3d nt = F*n0;
		
		// calculate the outer product of nt
		mat3ds N = dyad(nt);
		tens4ds NxN = dyad1s(N);
		
		// calculate strain energy derivative
		double Wl = m_ksi*pow(In_1, m_beta-1.0)*exp(m_alpha*pow(In_1, m_beta));
		
		// calculate the fiber stress
		s = N*(2.0*Wl/J);
        
		// calculate strain energy 2nd derivative
		double tmp = m_alpha*pow(In_1, m_beta);
		double Wll = m_ksi*pow(In_1, m_beta-2.0)*((tmp+1)*m_beta-1.0)*exp(tmp);
		
		// calculate the fiber tangent
		c = NxN*(4.0*Wll/J);
        
        // add the contribution from shear
        mat3ds B = pt.DevLeftCauchyGreen();
        c += dyad4s(N,B)*(m_mu/J);
	}
	else
	{
		c.zero();
	}
	
	// This is the final value of the elasticity tensor
    mat3dd I(1);
	tens4ds IxI = dyad1s(I);
	tens4ds I4  = dyad4s(I);
	c += ((I4+IxI/3.0)*s.tr() - dyad1s(I,s))*(2./3.) - (ddots(IxI, c)-IxI*(c.tr()/3.))/3.;
    
	return c;
}

//-----------------------------------------------------------------------------
//! Strain energy density
double FEFiberExponentialPowerUC::DevStrainEnergyDensity(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    FEFiberMaterialPoint& pf = *mp.ExtractData<FEFiberMaterialPoint>();
	
	// loop over all integration points
	mat3ds C = pt.DevRightCauchyGreen();
    mat3ds C2 = C*C;
	
	// fiber direction in global coordinate system
	vec3d n0 = GetFiberVector(mp);
	
	// Calculate In = n0*C*n0
	double In_1 = n0*(C*n0) - 1.0;
	
	// only take fibers in tension into consideration
	const double eps = 0;
	double sed = 0.0;
	if (In_1 > eps)
	{
		// calculate strain energy density
        if (m_alpha > 0)
            sed = m_ksi/(m_alpha*m_beta)*(exp(m_alpha*pow(In_1, m_beta))-1);
        else
            sed = m_ksi/m_beta*pow(In_1, m_beta);
		
        // add the contribution from shear
        sed += m_mu*(n0*(C2*n0)-2*In_1-1)/4.0;
	}
    
    return sed;
}

//-----------------------------------------------------------------------------
// FEFiberNH
//-----------------------------------------------------------------------------

// define the material parameters
BEGIN_PARAMETER_LIST(FEFiberNHUC, FEElasticFiberMaterialUC)
    ADD_PARAMETER2(m_mu, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "mu");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
mat3ds FEFiberNHUC::DevStress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	double J = pt.m_J;
	// distortional part of deformation gradient
	mat3d F = pt.m_F*pow(J,-1.0/3.0);
	
	// loop over all integration points
	const double eps = 0;
	mat3ds C = pt.DevRightCauchyGreen();
	mat3ds s;
	
	// fiber direction in global coordinate system
	vec3d n0 = GetFiberVector(mp);
	
	// Calculate In = n0*C*n0
	double In_1 = n0*(C*n0) - 1.0;
	
	// only take fibers in tension into consideration
	if (In_1 > eps)
	{
		// get the global spatial fiber direction in current configuration
		vec3d nt = F*n0;
		
		// calculate the outer product of nt
		mat3ds N = dyad(nt);
		
		// calculate the fiber stress
		s = N*(m_mu*In_1/J);
	}
	else
	{
		s.zero();
	}
	
	return s.dev();
}

//-----------------------------------------------------------------------------
tens4ds FEFiberNHUC::DevTangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	double J = pt.m_J;
	// distortional part of deformation gradient
	mat3d F = pt.m_F*pow(J,-1.0/3.0);
	
	// loop over all integration points
	const double eps = 0;
	mat3ds C = pt.DevRightCauchyGreen();
    mat3ds s;
	tens4ds c;
	
	// fiber direction in global coordinate system
	vec3d n0 = GetFiberVector(mp);
	
	// Calculate In = n0*C*n0
	double In_1 = n0*(C*n0) - 1.0;
	
	// only take fibers in tension into consideration
	if (In_1 > eps)
	{
		// get the global spatial fiber direction in current configuration
		vec3d nt = F*n0;
		
		// calculate the outer product of nt
		mat3ds N = dyad(nt);
		tens4ds NxN = dyad1s(N);
		
		// calculate the fiber stress
		s = N*(m_mu*In_1/J);
        
		// calculate the fiber tangent
		c = NxN*(2*m_mu/J);
	}
	else
	{
		c.zero();
	}
	
	// This is the final value of the elasticity tensor
    mat3dd I(1);
	tens4ds IxI = dyad1s(I);
	tens4ds I4  = dyad4s(I);
	c += ((I4+IxI/3.0)*s.tr() - dyad1s(I,s))*(2./3.)
	- (ddots(IxI, c)-IxI*(c.tr()/3.))/3.;

	return c;
}

//-----------------------------------------------------------------------------
//! Strain energy density
double FEFiberNHUC::DevStrainEnergyDensity(FEMaterialPoint& mp)
{
    double sed = 0.0;
    
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// loop over all integration points
	const double eps = 0;
	mat3ds C = pt.DevRightCauchyGreen();
	
	// fiber direction in global coordinate system
	vec3d n0 = GetFiberVector(mp);
	
	// Calculate In = n0*C*n0
	double In_1 = n0*(C*n0) - 1.0;
	
	// only take fibers in tension into consideration
	if (In_1 > eps)
        sed = 0.25*m_mu*In_1*In_1;
    
    return sed;
}
