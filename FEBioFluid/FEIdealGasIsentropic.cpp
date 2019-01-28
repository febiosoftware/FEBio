#include "FEIdealGasIsentropic.h"
#include "FECore/FEModel.h"
#include "FECore/FECoreKernel.h"
#include <FECore/fecore_error.h>

// define the material parameters
BEGIN_FECORE_CLASS(FEIdealGasIsentropic, FEFluid)
	ADD_PARAMETER(m_gamma, FE_RANGE_GREATER(0.0), "gamma");
	ADD_PARAMETER(m_M    , FE_RANGE_GREATER(0.0), "M"    );
END_FECORE_CLASS();

//============================================================================
// FEIdealGasIsentropic
//============================================================================

//-----------------------------------------------------------------------------
//! FEIdealGasIsentropic constructor

FEIdealGasIsentropic::FEIdealGasIsentropic(FEModel* pfem) : FEFluid(pfem)
{
    m_rhor = 0;
    m_k = 0;
    m_gamma = 0;
    m_M = 0;
}

//-----------------------------------------------------------------------------
//! initialization
bool FEIdealGasIsentropic::Init() 
{
    m_R  = GetFEModel()->GetGlobalConstant("R");
    m_Tr = GetFEModel()->GetGlobalConstant("T");
    m_pr = GetFEModel()->GetGlobalConstant("p");
    
    if (m_R <= 0) return fecore_error("A positive universal gas constant R must be defined in Globals section");
    if (m_Tr <= 0) return fecore_error("A positive ambient absolute temperature T must be defined in Globals section");
    if (m_pr <= 0) return fecore_error("A positive ambient absolute pressure p must be defined in Globals section");

    m_rhor = m_M*m_pr/(m_R*m_Tr);
    
    return true;
}

//-----------------------------------------------------------------------------
//! elastic pressure from dilatation
double FEIdealGasIsentropic::Pressure(const double e)
{
    double J = 1 + e;
    return m_pr*(pow(J, -m_gamma) - 1);
}

//-----------------------------------------------------------------------------
//! tangent of elastic pressure with respect to strain J
double FEIdealGasIsentropic::Tangent_Pressure_Strain(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    double J = fp.m_Jf;
    double dp = -m_gamma*m_pr*pow(J, -m_gamma-1);
    return dp;
}

//-----------------------------------------------------------------------------
//! 2nd tangent of elastic pressure with respect to strain J
double FEIdealGasIsentropic::Tangent_Pressure_Strain_Strain(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    double J = fp.m_Jf;
    double d2p = m_gamma*(m_gamma+1)*m_pr*pow(J, -m_gamma-2);
    return d2p;
}

//-----------------------------------------------------------------------------
//! evaluate temperature
double FEIdealGasIsentropic::Temperature(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    double J = fp.m_Jf;
    double T = m_Tr*pow(J, 1-m_gamma);
    return T;
}

//-----------------------------------------------------------------------------
//! calculate free energy density (per reference volume)
double FEIdealGasIsentropic::StrainEnergyDensity(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    double J = fp.m_Jf;
    double sed = m_pr*(J-1+(pow(J, 1-m_gamma)-1)/(m_gamma-1));
    return sed;
}

//-----------------------------------------------------------------------------
//! invert pressure-dilatation relation
double FEIdealGasIsentropic::Dilatation(const double p)
{
    double J = pow(p/m_pr+1, -1./m_gamma);
    return J - 1;
}

