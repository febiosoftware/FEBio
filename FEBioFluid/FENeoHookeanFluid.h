#pragma once
#include "FEElasticFluid.h"

//-----------------------------------------------------------------------------
// This class evaluates the pressure for a neo-Hookean fluid
// suitable for modeling nearly incompressible behavior

class FENeoHookeanFluid :	public FEElasticFluid
{
public:
    //! constructor
    FENeoHookeanFluid(FEModel* pfem);
    
    //! fluid pressure
    double Pressure(FEMaterialPoint& pt);
    
    //! tangent of fluid pressure with respect to strain J
    double Tangent_Pressure_Strain(FEMaterialPoint& mp);
    
public:
    double	m_k;		//!< bulk modulus
    
    // declare parameter list
    DECLARE_PARAMETER_LIST();
};
