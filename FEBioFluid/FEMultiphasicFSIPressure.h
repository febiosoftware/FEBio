//
//  FEFluidSolutesPressure.hpp
//  FEBioFluid
//
//  Created by Jay Shim on 12/10/20.
//  Copyright © 2020 febio.org. All rights reserved.
//

#pragma once
#include <FECore/FESurfaceLoad.h>
#include "FEMultiphasicFSI.h"

//-----------------------------------------------------------------------------
//! FEFluidResistanceBC is a fluid surface that has a normal
//! pressure proportional to the flow rate (resistance).
//!
class FEBIOFLUID_API FEMultiphasicFSIPressure : public FESurfaceLoad
{
public:
    //! constructor
    FEMultiphasicFSIPressure(FEModel* pfem);
    
    //! calculate traction stiffness (there is none)
    void StiffnessMatrix(FELinearSystem& LS) override {}
    
    //! calculate load vector
    void LoadVector(FEGlobalVector& R) override;
    
    //! set the dilatation
    void Update() override;
    
    
    //! initialize
    bool Init() override;
    
    //! activate
    void Activate() override;
    
    //! serialization
    void Serialize(DumpStream& ar) override;
    
private:
    double          m_p;       //!< prescribed fluid pressure
    
private:
    FEMultiphasicFSI*    m_pfs;   //!< pointer to fluid-solutes material
    
    int        m_dofEF;
    int        m_dofC;
    
    DECLARE_FECORE_CLASS();
};
