#include "FEIsotropicFourier.h"

//-----------------------------------------------------------------------------
// define the parameter list
BEGIN_PARAMETER_LIST(FEIsotropicFourier, FEMaterial)
	ADD_PARAMETER(m_k  , FE_PARAM_DOUBLE, "k");
	ADD_PARAMETER(m_c  , FE_PARAM_DOUBLE, "c");
	ADD_PARAMETER(m_rho, FE_PARAM_DOUBLE, "density");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Initialize isotropic Fourier material data
void FEIsotropicFourier::Init()
{
	if (m_k <= 0) throw MaterialError("Invalid value for k");
	if (m_c <= 0) throw MaterialError("Invalid value for c");
	if (m_rho <= 0) throw MaterialError("Invalid value for density");
}

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
