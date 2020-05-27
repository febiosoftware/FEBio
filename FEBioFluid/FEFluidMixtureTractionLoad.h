//
//  FEFluidMixtureTractionLoad.hpp
//  FEBioFluid
//
//  Created by Jay Shim on 3/2/20.
//  Copyright © 2020 febio.org. All rights reserved.
//

#pragma once
#include <FECore/FESurfaceLoad.h>
#include <FECore/FEModelParam.h>
#include "febiofluid_api.h"

//-----------------------------------------------------------------------------
//! FEFluidTractionLoad is a fluid surface that has a prescribed
//! viscous traction vector on it.
//!
class FEBIOFLUID_API FEFluidMixtureTractionLoad : public FESurfaceLoad
{
public:
    //! constructor
    FEFluidMixtureTractionLoad(FEModel* pfem);
    
    //! initialization
    bool Init() override;
    
    //! Set the surface to apply the load to
    void SetSurface(FESurface* ps) override;
    
    //! calculate traction stiffness (there is none)
    void StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp) override;
    
    //! calculate load vector
    void LoadVector(FEGlobalVector& R, const FETimeInfo& tp) override;
    
private:
    double            m_scale;    //!< magnitude of traction load
    FEParamVec3     m_TC;        //!< traction boundary cards
    
    DECLARE_FECORE_CLASS();
};
