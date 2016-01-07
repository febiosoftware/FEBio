#pragma once
#include "FEElasticFluid.h"

//-----------------------------------------------------------------------------
// This class evaluates the pressure for an ideal gas

class FEIdealGas :	public FEElasticFluid
	{
	public:
		//! constructor
		FEIdealGas(FEModel* pfem);
		
        //! fluid pressure
        double Pressure(FEMaterialPoint& pt);
        
        //! tangent of fluid pressure with respect to strain J
        double Tangent_Pressure_Strain(FEMaterialPoint& mp);
		
        //! 2nd derivative of fluid pressure with respect to strain J
        double Tangent_Pressure_Strain_Strain(FEMaterialPoint& mp);
        
        //! Initialization routine
        bool Init();
        
	public:
		double	m_M;		//!< molar mass
        double	m_Rgas;		//!< universal gas constant
        double	m_Tabs;		//!< absolute temperature
		
		// declare parameter list
		DECLARE_PARAMETER_LIST();
	};
