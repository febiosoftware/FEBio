//
//  FEIdealFluid.h
//  FEBioFluid
//
//  Created by Gerard Ateshian on 9/16/15.
//  Copyright (c) 2015 febio.org. All rights reserved.
//

#pragma once
#include "FEElasticFluid.h"

//-----------------------------------------------------------------------------
// This class evaluates the pressure for an ideal gas

class FEIdealFluid :	public FEElasticFluid
{
public:
    //! constructor
    FEIdealFluid(FEModel* pfem);
    
    //! fluid pressure
    double Pressure(FEMaterialPoint& pt);
    
    //! tangent of fluid pressure with respect to strain J
    double Tangent_Pressure_Strain(FEMaterialPoint& mp);
    
    //! Initialization routine
    void Init();
    
public:
    double	m_k;		//!< bulk modulus
    
    // declare parameter list
    DECLARE_PARAMETER_LIST();
};
