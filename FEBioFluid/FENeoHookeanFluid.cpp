#include "FENeoHookeanFluid.h"
#include "FEFluid.h"
#include "FECore/FEModel.h"

// define the material parameters
BEGIN_PARAMETER_LIST(FENeoHookeanFluid, FEElasticFluid)
ADD_PARAMETER2(m_k, FE_PARAM_DOUBLE, FE_RANGE_GREATER(0.0), "k");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Constructor.
FENeoHookeanFluid::FENeoHookeanFluid(FEModel* pfem) : FEElasticFluid(pfem)
{
    m_k = 0;
}

//-----------------------------------------------------------------------------
//! fluid pressure (relative to reference pressure when J=1)
double FENeoHookeanFluid::Pressure(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& vt = *mp.ExtractData<FEFluidMaterialPoint>();
    double J = vt.m_J;
    double p = -m_k*log(J)/J;
    return p;
}

//-----------------------------------------------------------------------------
//! tangent of fluid pressure with respect to strain J
double FENeoHookeanFluid::Tangent_Pressure_Strain(FEMaterialPoint &mp)
{
    FEFluidMaterialPoint& vt = *mp.ExtractData<FEFluidMaterialPoint>();
    double J = vt.m_J;
    double dpdJ = m_k*(log(J) - 1)/(J*J);
    return dpdJ;
}

//-----------------------------------------------------------------------------
//! 2nd derivative of fluid pressure with respect to strain J
double FENeoHookeanFluid::Tangent_Pressure_Strain_Strain(FEMaterialPoint &mp)
{
    FEFluidMaterialPoint& vt = *mp.ExtractData<FEFluidMaterialPoint>();
    double J = vt.m_J;
    double d2pdJ2 = m_k*(3 - 2*log(J))/pow(J, 3);
    return d2pdJ2;
}
