#include "FEElasticFiberMaterial.h"
#include "FEFiberMaterialPoint.h"

//-----------------------------------------------------------------------------
// FEFiberExponentialPower
//-----------------------------------------------------------------------------

// define the material parameters
BEGIN_PARAMETER_LIST(FEFiberExponentialPower, FEElasticFiberMaterial)
	ADD_PARAMETER(m_alpha, FE_PARAM_DOUBLE, "alpha");
	ADD_PARAMETER(m_beta, FE_PARAM_DOUBLE, "beta");
	ADD_PARAMETER(m_ksi , FE_PARAM_DOUBLE, "ksi" );
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
void FEFiberExponentialPower::Init()
{
    FEMaterial::Init();
	if (m_ksi < 0) throw MaterialError("ksi must be positive.");
	if (m_beta < 2) throw MaterialError("beta must be >= 2.");
	if (m_alpha < 0) throw MaterialError("alpha must be >= 0.");
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
	}
	else
	{
		c.zero();
	}
	
	return c;
}

//-----------------------------------------------------------------------------
// FEFiberNH
//-----------------------------------------------------------------------------

// define the material parameters
BEGIN_PARAMETER_LIST(FEFiberNH, FEElasticFiberMaterial)
	ADD_PARAMETER(m_mu, FE_PARAM_DOUBLE, "mu");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
void FEFiberNH::Init()
{
    FEMaterial::Init();
	if (m_mu < 0) throw MaterialError("mu must be positive.");
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
