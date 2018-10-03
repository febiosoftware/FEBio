//
//  FEActiveConstantSupply.h
//  FEBioMix
//
//  Created by Gerard Ateshian on 12/9/14.
//  Copyright (c) 2014 febio.org. All rights reserved.
//

#ifndef __FEBioMix__FEActiveConstantSupply__
#define __FEBioMix__FEActiveConstantSupply__

#include "FEActiveMomentumSupply.h"

//-----------------------------------------------------------------------------
// This class implements a constant active momentum supply

class FEActiveConstantSupply :	public FEActiveMomentumSupply
{
public:
    //! constructor
    FEActiveConstantSupply(FEModel* pfem);
    
    //! active momentum supply
    vec3d ActiveSupply(FEMaterialPoint& pt) override;
    
    //! Tangent of active momentum supply
    vec3d Tangent_ActiveSupply_Strain(FEMaterialPoint& mp) override;
        
public:
    double	m_asupp;			//!< active momentum supply
    
    // declare parameter list
    DECLARE_FECORE_CLASS();
};

#endif /* defined(__FEBioMix__FEActiveConstantSupply__) */
