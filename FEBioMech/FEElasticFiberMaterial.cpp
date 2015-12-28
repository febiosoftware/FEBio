#include "FEElasticFiberMaterial.h"
#include "FEFiberMaterialPoint.h"

//-----------------------------------------------------------------------------
// FEFiberExponentialPower
//-----------------------------------------------------------------------------

// define the material parameters
BEGIN_PARAMETER_LIST(FEFiberExponentialPower, FEElasticFiberMaterial)
	ADD_PARAMETER2(m_alpha, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "alpha");
	ADD_PARAMETER2(m_beta , FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(2.0), "beta" );
	ADD_PARAMETER(m_ksi  , FE_PARAM_DOUBLE, "ksi"  );
    ADD_PARAMETER2(m_mu   , FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "mu"   );
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
bool FEFiberExponentialPower::Init()
{
	if ((4*m_ksi + 2*m_mu) < 0) return MaterialError("4*ksi+2*mu must be positive.");
    return FEElasticFiberMaterial::Init();
}

//-----------------------------------------------------------------------------
void FEElasticFiberMaterial::SetFiberDirection(FEMaterialPoint& mp, const vec3d n0)
{
    FEFiberMaterialPoint& pf = *mp.ExtractData<FEFiberMaterialPoint>();
    pf.m_n0 = n0;
}

//-----------------------------------------------------------------------------
mat3ds FEFiberExponentialPower::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    FEFiberMaterialPoint& pf = *mp.ExtractData<FEFiberMaterialPoint>();
	
	// deformation gradient
	mat3d &F = pt.m_F;
	double J = pt.m_J;
	
	// loop over all integration points
	vec3d n0, nt;
	double In_1, Wl;
	const double eps = 0;
	mat3ds C = pt.RightCauchyGreen();
	mat3ds s;
	
	// fiber direction in global coordinate system
	n0 = pf.m_n0;
	
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
        
        // add the contribution from shear
        mat3ds BmI = pt.LeftCauchyGreen() - mat3dd(1);
        s += (N*BmI + BmI*N)*(m_mu/2/J);
	}
	else
	{
		s.zero();
	}
	
	return s;
}

//-----------------------------------------------------------------------------
tens4ds FEFiberExponentialPower::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    FEFiberMaterialPoint& pf = *mp.ExtractData<FEFiberMaterialPoint>();
	
	// deformation gradient
	mat3d &F = pt.m_F;
	double J = pt.m_J;
	
	vec3d n0, nt;
	double In_1, Wll;
	const double eps = 0;
	mat3ds C = pt.RightCauchyGreen();
	tens4ds c;
	
	// fiber direction in global coordinate system
	n0 = pf.m_n0;
	
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
        
        // add the contribution from shear
        mat3ds B = pt.LeftCauchyGreen();
        c += dyad4s(N,B)*(m_mu/J);
	}
	else
	{
		c.zero();
	}
	
	return c;
}

//-----------------------------------------------------------------------------
//! Strain energy density
double FEFiberExponentialPower::StrainEnergyDensity(FEMaterialPoint& mp)
{
    double sed = 0.0;
    
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    FEFiberMaterialPoint& pf = *mp.ExtractData<FEFiberMaterialPoint>();
	
	// loop over all integration points
	vec3d n0, nt;
	double In_1;
	const double eps = 0;
	mat3ds C = pt.RightCauchyGreen();
    mat3ds C2 = C*C;
	
	// fiber direction in global coordinate system
	n0 = pf.m_n0;
	
	// Calculate In = n0*C*n0
	In_1 = n0*(C*n0) - 1.0;
	
	// only take fibers in tension into consideration
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
BEGIN_PARAMETER_LIST(FEFiberNH, FEElasticFiberMaterial)
	ADD_PARAMETER(m_mu, FE_PARAM_DOUBLE, "mu");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
bool FEFiberNH::Init()
{
	if (m_mu < 0) return MaterialError("mu must be positive.");
    return FEElasticFiberMaterial::Init();
}

//-----------------------------------------------------------------------------
mat3ds FEFiberNH::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    FEFiberMaterialPoint& pf = *mp.ExtractData<FEFiberMaterialPoint>();
	
	// deformation gradient
	mat3d &F = pt.m_F;
	double J = pt.m_J;
	
	// loop over all integration points
	vec3d n0, nt;
	double In_1;
	const double eps = 0;
	mat3ds C = pt.RightCauchyGreen();
	mat3ds s;
	
	// fiber direction in global coordinate system
	n0 = pf.m_n0;
	
	// Calculate In = n0*C*n0
	In_1 = n0*(C*n0) - 1.0;
	
	// only take fibers in tension into consideration
	if (In_1 > eps)
	{
		// get the global spatial fiber direction in current configuration
		nt = F*n0;
		
		// calculate the outer product of nt
		mat3ds N = dyad(nt);
		
		// calculate the fiber stress
		s = N*(m_mu*In_1/J);
	}
	else
	{
		s.zero();
	}
	
	return s;
}

//-----------------------------------------------------------------------------
tens4ds FEFiberNH::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    FEFiberMaterialPoint& pf = *mp.ExtractData<FEFiberMaterialPoint>();
	
	// deformation gradient
	mat3d &F = pt.m_F;
	double J = pt.m_J;
	
	// loop over all integration points
	vec3d n0, nt;
	double In_1;
	const double eps = 0;
	mat3ds C = pt.RightCauchyGreen();
	tens4ds c;
	
	// fiber direction in global coordinate system
	n0 = pf.m_n0;
	
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
		
		// calculate the fiber tangent
		c = NxN*(2*m_mu/J);
	}
	else
	{
		c.zero();
	}
	
	return c;
}

//-----------------------------------------------------------------------------
//! Strain energy density
double FEFiberNH::StrainEnergyDensity(FEMaterialPoint& mp)
{
    double sed = 0.0;
    
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    FEFiberMaterialPoint& pf = *mp.ExtractData<FEFiberMaterialPoint>();
	
	// loop over all integration points
	vec3d n0, nt;
	double In_1;
	const double eps = 0;
	mat3ds C = pt.RightCauchyGreen();
	
	// fiber direction in global coordinate system
	n0 = pf.m_n0;
	
	// Calculate In = n0*C*n0
	In_1 = n0*(C*n0) - 1.0;
	
	// only take fibers in tension into consideration
	if (In_1 > eps)
        sed = 0.25*m_mu*In_1*In_1;
    
    return sed;
}

//-----------------------------------------------------------------------------
// FEFiberPowerLinear
//-----------------------------------------------------------------------------

// define the material parameters
BEGIN_PARAMETER_LIST(FEFiberPowerLinear, FEElasticFiberMaterial)
    ADD_PARAMETER(m_E    , FE_PARAM_DOUBLE, "E"    );
    ADD_PARAMETER(m_beta , FE_PARAM_DOUBLE, "beta" );
    ADD_PARAMETER(m_lam0 , FE_PARAM_DOUBLE, "lam0" );
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
bool FEFiberPowerLinear::Init()
{
    if (FEElasticFiberMaterial::Init() == false) return false;
    if (m_E < 0) return MaterialError("E must be positive.");
    if (m_beta < 2) return MaterialError("beta must be >= 2.");
    if (m_lam0 <= 1) return MaterialError("lam0 must be >1.");
    
    // initialize material constants
    m_I0 = m_lam0*m_lam0;
    m_ksi = m_E/4/(m_beta-1)*pow(m_I0, -3./2.)*pow(m_I0-1, 2-m_beta);
    m_b = m_ksi*pow(m_I0-1, m_beta-1) + m_E/2/sqrt(m_I0);

	return true;
}

//-----------------------------------------------------------------------------
mat3ds FEFiberPowerLinear::Stress(FEMaterialPoint& mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    FEFiberMaterialPoint& pf = *mp.ExtractData<FEFiberMaterialPoint>();
    
    // deformation gradient
    mat3d &F = pt.m_F;
    double J = pt.m_J;
    
    // loop over all integration points
    vec3d n0, nt;
    double In, sn;
    const double eps = 0;
    mat3ds C = pt.RightCauchyGreen();
    mat3ds s;
    
    // fiber direction in global coordinate system
    n0 = pf.m_n0;
    
    // Calculate In
    In = n0*(C*n0);
    
    // only take fibers in tension into consideration
    if (In - 1 > eps)
    {
        // get the global spatial fiber direction in current configuration
        nt = F*n0/sqrt(In);
        
        // calculate the outer product of nt
        mat3ds N = dyad(nt);
        
        // calculate the fiber stress magnitude
        sn = (In < m_I0) ?
        2*In*m_ksi*pow(In-1, m_beta-1) :
        2*m_b*In - m_E*sqrt(In);
        
        // calculate the fiber stress
        s = N*(sn/J);
    }
    else
    {
        s.zero();
    }
    
    return s;
}

//-----------------------------------------------------------------------------
tens4ds FEFiberPowerLinear::Tangent(FEMaterialPoint& mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    FEFiberMaterialPoint& pf = *mp.ExtractData<FEFiberMaterialPoint>();
    
    // deformation gradient
    mat3d &F = pt.m_F;
    double J = pt.m_J;
    
    vec3d n0, nt;
    double In, cn;
    const double eps = 0;
    mat3ds C = pt.RightCauchyGreen();
    tens4ds c;
    
    // fiber direction in global coordinate system
    n0 = pf.m_n0;
    
    // Calculate In
    In = n0*(C*n0);
    
    // only take fibers in tension into consideration
    if (In - 1 > eps)
    {
        // get the global spatial fiber direction in current configuration
        nt = F*n0/sqrt(In);
        
        // calculate the outer product of nt
        mat3ds N = dyad(nt);
        tens4ds NxN = dyad1s(N);
        
        // calculate modulus
        cn = (In < m_I0) ?
        4*In*In*m_ksi*(m_beta-1)*pow(In-1, m_beta-2) :
        m_E*sqrt(In);
        
        // calculate the fiber tangent
        c = NxN*(cn/J);
    }
    else
    {
        c.zero();
    }
    
    return c;
}

//-----------------------------------------------------------------------------
//! Strain energy density
double FEFiberPowerLinear::StrainEnergyDensity(FEMaterialPoint& mp)
{
    double sed = 0.0;
    
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    FEFiberMaterialPoint& pf = *mp.ExtractData<FEFiberMaterialPoint>();
    
    // loop over all integration points
    vec3d n0, nt;
    double In;
    const double eps = 0;
    mat3ds C = pt.RightCauchyGreen();
    
    // fiber direction in global coordinate system
    n0 = pf.m_n0;
    
    // Calculate In = n0*C*n0
    In = n0*(C*n0);
    
    // only take fibers in tension into consideration
    if (In - 1 > eps)
    {
        // calculate strain energy density
        sed = (In < m_I0) ?
        m_ksi/m_beta*pow(In-1, m_beta) :
        m_b*(In-m_I0) - m_E*(sqrt(In)-sqrt(m_I0)) + m_ksi/m_beta*pow(m_I0-1, m_beta);
    }
    
    return sed;
}

