#include "FEFiberPowLinear.h"

// define the material parameters
BEGIN_PARAMETER_LIST(FEFiberPowLinear, FEElasticFiberMaterial)
    ADD_PARAMETER2(m_E    , FE_PARAM_DOUBLE, FE_RANGE_GREATER(0.0), "E"    );
    ADD_PARAMETER2(m_lam0 , FE_PARAM_DOUBLE, FE_RANGE_GREATER(1.0), "lam0" );
    ADD_PARAMETER2(m_beta , FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(2.0), "beta" );
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
// FEFiberPowLinear
//-----------------------------------------------------------------------------

FEFiberPowLinear::FEFiberPowLinear(FEModel* pfem) : FEElasticFiberMaterial(pfem)
{

}

//-----------------------------------------------------------------------------
mat3ds FEFiberPowLinear::Stress(FEMaterialPoint& mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // initialize material constants
    double I0 = m_lam0*m_lam0;
    double ksi = m_E/4/(m_beta-1)*pow(I0, -3./2.)*pow(I0-1, 2-m_beta);
    double b = ksi*pow(I0-1, m_beta-1) + m_E/2/sqrt(I0);
    
    // deformation gradient
    mat3d &F = pt.m_F;
    double J = pt.m_J;
    
    // loop over all integration points
    const double eps = 0;
    mat3ds C = pt.RightCauchyGreen();
    mat3ds s;
    
    // evaluate fiber direction in global coordinate system
    vec3d n0 = GetFiberVector(mp);
    
    // Calculate In
    double In = n0*(C*n0);
    
    // only take fibers in tension into consideration
    if (In - 1 > eps)
    {
        // get the global spatial fiber direction in current configuration
        vec3d nt = F*n0/sqrt(In);
        
        // calculate the outer product of nt
        mat3ds N = dyad(nt);
        
        // calculate the fiber stress magnitude
        double sn = (In < I0) ?
        2*In*ksi*pow(In-1, m_beta-1) :
        2*b*In - m_E*sqrt(In);
        
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
tens4ds FEFiberPowLinear::Tangent(FEMaterialPoint& mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // initialize material constants
    double I0 = m_lam0*m_lam0;
    double ksi = m_E/4/(m_beta-1)*pow(I0, -3./2.)*pow(I0-1, 2-m_beta);
    
    // deformation gradient
    mat3d &F = pt.m_F;
    double J = pt.m_J;
    
    // loop over all integration points
    const double eps = 0;
    mat3ds C = pt.RightCauchyGreen();
    tens4ds c;
    
    // evaluate fiber direction in global coordinate system
    vec3d n0 = GetFiberVector(mp);
    
    // Calculate In
    double In = n0*(C*n0);
    
    // only take fibers in tension into consideration
    if (In - 1 > eps)
    {
        // get the global spatial fiber direction in current configuration
        vec3d nt = F*n0/sqrt(In);
        
        // calculate the outer product of nt
        mat3ds N = dyad(nt);
        tens4ds NxN = dyad1s(N);
        
        // calculate modulus
        double cn = (In < I0) ?
        4*In*In*ksi*(m_beta-1)*pow(In-1, m_beta-2) :
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
double FEFiberPowLinear::StrainEnergyDensity(FEMaterialPoint& mp)
{
    double sed = 0.0;
    
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // initialize material constants
    double I0 = m_lam0*m_lam0;
    double ksi = m_E/4/(m_beta-1)*pow(I0, -3./2.)*pow(I0-1, 2-m_beta);
    double b = ksi*pow(I0-1, m_beta-1) + m_E/2/sqrt(I0);
    
    // loop over all integration points
    const double eps = 0;
    mat3ds C = pt.RightCauchyGreen();
    
    // evaluate fiber direction in global coordinate system
    vec3d n0 = GetFiberVector(mp);
    
    // Calculate In = n0*C*n0
    double In = n0*(C*n0);
    
    // only take fibers in tension into consideration
    if (In - 1 > eps)
    {
        // calculate strain energy density
        sed = (In < I0) ?
        ksi/m_beta*pow(In-1, m_beta) :
        b*(In-I0) - m_E*(sqrt(In)-sqrt(I0)) + ksi/m_beta*pow(I0-1, m_beta);
    }
    
    return sed;
}
