//
//  FEBackFlowBiphasicStabilization.hpp
//  FEBioFluid
//
//  Created by Jay Shim on 2/6/20.
//  Copyright © 2020 febio.org. All rights reserved.
//

#pragma once
#include <FECore/FESurfaceLoad.h>
#include "febiofluid_api.h"

//-----------------------------------------------------------------------------
//! Backflow stabilization prescribes a normal traction that opposes
//! backflow on a boundary surface.
class FEBIOFLUID_API FEBackFlowBiphasicStabilization : public FESurfaceLoad
{
public:
    //! constructor
    FEBackFlowBiphasicStabilization(FEModel* pfem);
    
    //! calculate pressure stiffness
    void StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp) override;
    
    //! calculate residual
    void LoadVector(FEGlobalVector& R, const FETimeInfo& tp) override;
    
    //! serialize data
    void Serialize(DumpStream& ar) override;
    
    //! initialization
    bool Init() override;
    
protected:
    vec3d FluidVelocity(FESurfaceMaterialPoint& mp, double alpha);
    
protected:
    double            m_beta;     //!< backflow stabilization coefficient
    double          m_rho;      //!< fluid density
    
    // degrees of freedom
    FEDofList    m_dofU;
    FEDofList    m_dofW;
    
    DECLARE_FECORE_CLASS();
};
