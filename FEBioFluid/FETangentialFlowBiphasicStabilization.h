//
//  FETangentialFlowBiphasicStabilization.hpp
//  FEBioFluid
//
//  Created by Jay Shim on 2/6/20.
//  Copyright © 2020 febio.org. All rights reserved.
//

#pragma once
#include <FECore/FESurfaceLoad.h>
#include "febiofluid_api.h"

//-----------------------------------------------------------------------------
//! Tangential flow stabilization prescribes a shear traction that opposes
//! tangential fluid velocity on a boundary surface, in the presence of normal
//! flow.  This can help stabilize inflow/outflow conditions.
class FEBIOFLUID_API FETangentialFlowBiphasicStabilization : public FESurfaceLoad
{
public:
    //! constructor
    FETangentialFlowBiphasicStabilization(FEModel* pfem);
    
    //! Initialization
    bool Init() override;
    
    //! Set the surface to apply the load to
    void SetSurface(FESurface* ps) override;
    
    //! calculate pressure stiffness
    void StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp) override;
    
    //! calculate load vector
    void LoadVector(FEGlobalVector& R, const FETimeInfo& tp) override;
    
    //! serialize data
    void Serialize(DumpStream& ar) override;
    
protected:
    vec3d FluidVelocity(FESurfaceMaterialPoint& mp, double alpha);
    
protected:
    double            m_beta;     //!< damping coefficient
    double          m_rho;      //!< fluid density
    
    // degrees of freedom
    FEDofList    m_dofU;
    FEDofList    m_dofW;
    
    DECLARE_FECORE_CLASS();
};
