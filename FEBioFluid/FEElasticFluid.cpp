#include "FEElasticFluid.h"
#include "FEFluid.h"

double FEElasticFluid::BulkModulus(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& pt = *mp.ExtractData<FEFluidMaterialPoint>();

    double k = -pt.m_J*Tangent_Pressure_Strain(mp);

    return k;
}
