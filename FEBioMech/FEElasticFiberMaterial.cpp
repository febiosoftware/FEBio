#include "FEElasticFiberMaterial.h"
#include "FEFiberMaterialPoint.h"

BEGIN_PARAMETER_LIST(FEElasticFiberMaterial, FEElasticMaterial)
	ADD_PARAMETER(m_thd, FE_PARAM_DOUBLE, "theta");
	ADD_PARAMETER(m_phd, FE_PARAM_DOUBLE, "phi");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEElasticFiberMaterial::FEElasticFiberMaterial(FEModel* pfem) : FEElasticMaterial(pfem) 
{
	m_thd = 0.0;
	m_phd = 90.0;
}

//-----------------------------------------------------------------------------
FEMaterialPoint* FEElasticFiberMaterial::CreateMaterialPointData()
{
	FEFiberMaterialPoint* fp = new FEFiberMaterialPoint(FEElasticMaterial::CreateMaterialPointData());

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
vec3d FEElasticFiberMaterial::GetFiberVector(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	FEFiberMaterialPoint& fp = *mp.ExtractData<FEFiberMaterialPoint>();

	return pt.m_Q*fp.m_n0;
}

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
FEFiberExponentialPower::FEFiberExponentialPower(FEModel* pfem) : FEElasticFiberMaterial(pfem) 
{
	m_alpha = 0; 
	m_beta = 2; 
	m_ksi = 0; 
	m_mu = 0;
}

//-----------------------------------------------------------------------------
bool FEFiberExponentialPower::Validate()
{
	if ((4*m_ksi + 2*m_mu) < 0) return MaterialError("4*ksi+2*mu must be positive.");
    return FEElasticFiberMaterial::Validate();
}

//-----------------------------------------------------------------------------
mat3ds FEFiberExponentialPower::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    	
	// deformation gradient
	mat3d &F = pt.m_F;
	double J = pt.m_J;
	
	// loop over all integration points
	mat3ds C = pt.RightCauchyGreen();
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
		double Ib = pow(In_1, m_beta - 1.0);
		double Wl = m_ksi*Ib*exp(m_alpha*(Ib*In_1));
//		Wl = m_ksi*pow(In_1, m_beta-1.0)*exp(m_alpha*pow(In_1, m_beta));
		
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
	
	// deformation gradient
	mat3d &F = pt.m_F;
	double J = pt.m_J;
	
	mat3ds C = pt.RightCauchyGreen();
	tens4ds c;
	
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
		tens4ds NxN = dyad1s(N);
		
		// calculate strain energy 2nd derivative
		double tmp = m_alpha*pow(In_1, m_beta);
		double Wll = m_ksi*pow(In_1, m_beta-2.0)*((tmp+1)*m_beta-1.0)*exp(tmp);
		
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
	
	// loop over all integration points
	mat3ds C = pt.RightCauchyGreen();
    mat3ds C2 = C*C;
	
	// fiber direction in global coordinate system
	vec3d n0 = GetFiberVector(mp);
	
	// Calculate In = n0*C*n0
	double In_1 = n0*(C*n0) - 1.0;
	
	// only take fibers in tension into consideration
	const double eps = 0;
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
	ADD_PARAMETER2(m_mu, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "mu");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEFiberNH::FEFiberNH(FEModel* pfem) : FEElasticFiberMaterial(pfem) 
{ 
	m_mu = 0; 
}

//-----------------------------------------------------------------------------
mat3ds FEFiberNH::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// deformation gradient
	mat3d &F = pt.m_F;
	double J = pt.m_J;
	
	// loop over all integration points
	mat3ds C = pt.RightCauchyGreen();
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
	
	// deformation gradient
	mat3d &F = pt.m_F;
	double J = pt.m_J;
	
	// loop over all integration points
	mat3ds C = pt.RightCauchyGreen();
	tens4ds c;
	
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
	
	// loop over all integration points
	mat3ds C = pt.RightCauchyGreen();
	
	// fiber direction in global coordinate system
	vec3d n0 = GetFiberVector(mp);
	
	// Calculate In = n0*C*n0
	double In_1 = n0*(C*n0) - 1.0;
	
	// only take fibers in tension into consideration
	const double eps = 0;
	if (In_1 > eps)
        sed = 0.25*m_mu*In_1*In_1;
    
    return sed;
}

//-----------------------------------------------------------------------------
// FEFiberPowerLinear
//-----------------------------------------------------------------------------

// define the material parameters
BEGIN_PARAMETER_LIST(FEFiberPowerLinear, FEElasticFiberMaterial)
    ADD_PARAMETER2(m_E    , FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "E"    );
    ADD_PARAMETER2(m_beta , FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(2.0), "beta" );
    ADD_PARAMETER2(m_lam0 , FE_PARAM_DOUBLE, FE_RANGE_GREATER(1.0), "lam0" );
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEFiberPowerLinear::FEFiberPowerLinear(FEModel* pfem) : FEElasticFiberMaterial(pfem) 
{
	m_E = 0; 
	m_lam0 = 1; 
	m_beta = 3;
}

//-----------------------------------------------------------------------------
bool FEFiberPowerLinear::Validate()
{
    if (FEElasticFiberMaterial::Validate() == false) return false;
    
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
    
    // deformation gradient
    mat3d &F = pt.m_F;
    double J = pt.m_J;
    
    // loop over all integration points
    mat3ds C = pt.RightCauchyGreen();
    mat3ds s;
    
    // fiber direction in global coordinate system
	vec3d n0 = GetFiberVector(mp);
    
    // Calculate In
    double In = n0*(C*n0);
    
    // only take fibers in tension into consideration
	const double eps = 0;
	if (In - 1 > eps)
    {
        // get the global spatial fiber direction in current configuration
        vec3d nt = F*n0/sqrt(In);
        
        // calculate the outer product of nt
        mat3ds N = dyad(nt);
        
        // calculate the fiber stress magnitude
        double sn = (In < m_I0) ?
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
    
    // deformation gradient
    mat3d &F = pt.m_F;
    double J = pt.m_J;
    
    mat3ds C = pt.RightCauchyGreen();
    tens4ds c;
    
    // fiber direction in global coordinate system
	vec3d n0 = GetFiberVector(mp);
    
    // Calculate In
    double In = n0*(C*n0);
    
    // only take fibers in tension into consideration
	const double eps = 0;
	if (In - 1 > eps)
    {
        // get the global spatial fiber direction in current configuration
        vec3d nt = F*n0/sqrt(In);
        
        // calculate the outer product of nt
        mat3ds N = dyad(nt);
        tens4ds NxN = dyad1s(N);
        
        // calculate modulus
        double cn = (In < m_I0) ?
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
    
    // loop over all integration points
    mat3ds C = pt.RightCauchyGreen();
    
    // fiber direction in global coordinate system
	vec3d n0 = GetFiberVector(mp);
    
    // Calculate In = n0*C*n0
    double In = n0*(C*n0);
    
    // only take fibers in tension into consideration
	const double eps = 0;
	if (In - 1 > eps)
    {
        // calculate strain energy density
        sed = (In < m_I0) ?
        m_ksi/m_beta*pow(In-1, m_beta) :
        m_b*(In-m_I0) - m_E*(sqrt(In)-sqrt(m_I0)) + m_ksi/m_beta*pow(m_I0-1, m_beta);
    }
    
    return sed;
}
