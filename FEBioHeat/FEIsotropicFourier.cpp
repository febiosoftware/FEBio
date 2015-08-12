#include "FEIsotropicFourier.h"

//-----------------------------------------------------------------------------
// define the parameter list
BEGIN_PARAMETER_LIST(FEIsotropicFourier, FEMaterial)
	ADD_PARAMETER2(m_k  , FE_PARAM_DOUBLE, FE_RANGE_GREATER(0.0), "k");
	ADD_PARAMETER2(m_c  , FE_PARAM_DOUBLE, FE_RANGE_GREATER(0.0), "c");
	ADD_PARAMETER2(m_rho, FE_PARAM_DOUBLE, FE_RANGE_GREATER(0.0), "density");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
mat3ds FEIsotropicFourier::Conductivity(FEMaterialPoint& mp)
{
	return mat3dd(m_k);
}

//-----------------------------------------------------------------------------
void FEIsotropicFourier::Serialize(DumpFile& ar)
{
	FEHeatTransferMaterial::Serialize(ar);
}
