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
//! FluidV material point class.
//
class FEBIOFLUID_API FEFluidVMaterialPoint : public FEMaterialPoint
{
public:
    //! constructor
    FEFluidVMaterialPoint(FEMaterialPoint* pt);
    
    //! create a shallow copy
    FEMaterialPoint* Copy();
    
    //! data serialization
    void Serialize(DumpStream& ar);
    
    //! Data initialization
    void Init();
    
public:
    // FSI material data
    double          m_Jft;      //!< Jf at current time
    double          m_Jfp;      //!< Jf at previous time
    double          m_fJ;
    double          m_kJJ;
    double          m_dpf;      //!< derivative of fluid pressure w.r.t. Jf
    vector<vec3d>   m_kJv;
};

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
