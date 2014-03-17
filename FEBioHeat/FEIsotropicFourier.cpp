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
void FEIsotropicFourier::Conductivity(double D[3][3])
{
	D[0][0] = D[1][1] = D[2][2] = m_k;
	D[0][1] = D[1][0] = 0;
	D[0][2] = D[2][0] = 0;
	D[1][2] = D[2][1] = 0;
}

//-----------------------------------------------------------------------------
void FEIsotropicFourier::Serialize(DumpFile& ar)
{
	FEHeatTransferMaterial::Serialize(ar);
}
