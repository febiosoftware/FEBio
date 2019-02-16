//
//  FEFluidP.hpp
//  FEBioFluid
//
//  Created by Gerard Ateshian on 2/16/19.
//  Copyright © 2019 febio.org. All rights reserved.
//

#ifndef FEFluidP_hpp
#define FEFluidP_hpp

#include "FEFluid.h"

//-----------------------------------------------------------------------------
//! Fluids material is identical to fluid material

//-----------------------------------------------------------------------------
//! Base class for fluidP materials.

class FEBIOFLUID_API FEFluidP : public FEMaterial
{
public:
    FEFluidP(FEModel* pfem);
    
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
