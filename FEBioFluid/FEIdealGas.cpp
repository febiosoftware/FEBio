#include "FEIdealGas.h"
#include "FEFluid.h"
#include "FECore/FEModel.h"

// define the material parameters
BEGIN_PARAMETER_LIST(FEIdealGas, FEElasticFluid)
	ADD_PARAMETER2(m_M, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "molar_mass");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! Constructor. 
FEIdealGas::FEIdealGas(FEModel* pfem) : FEElasticFluid(pfem)
{
	m_M = 1;
}

//-----------------------------------------------------------------------------
// FEIdealGas
bool FEIdealGas::Init()
{
    m_Rgas = GetFEModel()->GetGlobalConstant("R");
    m_Tabs = GetFEModel()->GetGlobalConstant("T");
    
    return FEElasticFluid::Init();
}

//-----------------------------------------------------------------------------
//! fluid pressure (relative to reference pressure when J=1)
double FEIdealGas::Pressure(FEMaterialPoint& mp)
{
    FEFluid* m_pMat = dynamic_cast<FEFluid*>(GetParent());
    assert(m_pMat);
    double rho = m_pMat->Density(mp);
    double rhor = m_pMat->ReferentialDensity();
	return m_Rgas*m_Tabs*(rho - rhor)/m_M;
}

//-----------------------------------------------------------------------------
//! tangent of fluid pressure with respect to strain J
double FEIdealGas::Tangent_Pressure_Strain(FEMaterialPoint &mp)
{
    FEFluidMaterialPoint& vt = *mp.ExtractData<FEFluidMaterialPoint>();
    FEFluid* m_pMat = dynamic_cast<FEFluid*>(GetParent());
    assert(m_pMat);
    double rho = m_pMat->Density(mp);
    return -m_Rgas*m_Tabs*rho/m_M/vt.m_J;
}

//-----------------------------------------------------------------------------
//! 2nd derivative of fluid pressure with respect to strain J
double FEIdealGas::Tangent_Pressure_Strain_Strain(FEMaterialPoint &mp)
{
    FEFluidMaterialPoint& vt = *mp.ExtractData<FEFluidMaterialPoint>();
    FEFluid* m_pMat = dynamic_cast<FEFluid*>(GetParent());
    assert(m_pMat);
    double rho = m_pMat->Density(mp);
    return 2*m_Rgas*m_Tabs*rho/m_M/(vt.m_J*vt.m_J);
}
