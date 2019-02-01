//
//  FEFluidV.hpp
//  FEBioFluid
//
//  Created by Gerard Ateshian on 2/1/19.
//  Copyright © 2019 febio.org. All rights reserved.
//

#ifndef FEFluidV_hpp
#define FEFluidV_hpp

#include "FEFluid.h"

//-----------------------------------------------------------------------------
//! Fluids material is identical to fluid material

//-----------------------------------------------------------------------------
//! Base class for fluidV materials.

class FEBIOFLUID_API FEFluidV : public FEMaterial
{
public:
    FEFluidV(FEModel* pfem);
    
    // returns a pointer to a new material point object
    FEMaterialPoint* CreateMaterialPointData() override;
    
public:
    //! initialization
    bool Init() override;
    
public:
    FEFluid* Fluid() { return m_pFluid; }
    
private: // material properties
    FEFluid*                    m_pFluid;    //!< pointer to fluid material

    // declare parameter list
    DECLARE_FECORE_CLASS();
};

#endif /* FEFluidV_hpp */
