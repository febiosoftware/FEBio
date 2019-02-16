//
//  FEFluidPResistanceBC.hpp
//  FEBioFluid
//
//  Created by Gerard Ateshian on 2/16/19.
//  Copyright © 2019 febio.org. All rights reserved.
//

#ifndef FEFluidPResistanceBC_hpp
#define FEFluidPResistanceBC_hpp

#include <FECore/FESurfaceLoad.h>
#include <FECore/FESurfaceMap.h>
#include "FEFluid.h"

//-----------------------------------------------------------------------------
//! FEFluidPResistanceBC is a fluid surface that has a normal
//! pressure proportional to the flow rate (resistance).
//!
class FEBIOFLUID_API FEFluidPResistanceBC : public FESurfaceLoad
{
public:
    //! constructor
    FEFluidPResistanceBC(FEModel* pfem);
    
    //! Set the surface to apply the load to
    void SetSurface(FESurface* ps) override;
    
    //! calculate traction stiffness (there is none)
    void StiffnessMatrix(const FETimeInfo& tp, FESolver* psolver) override {}
    
    //! calculate residual
    void Residual(const FETimeInfo& tp, FEGlobalVector& R) override { m_alpha = tp.alpha; m_alphaf = tp.alphaf; }
    
    //! set the dilatation
    void Update() override;
    
    //! evaluate flow rate
    double FlowRate();
    
    //! initialize
    bool Init() override;
    
    //! activate
    void Activate() override;
    
private:
    double          m_R;        //!< flow resistance
    double          m_alpha;
    double          m_alphaf;
    double          m_p0;       //!< fluid pressure offset
    
    int        m_dofWX, m_dofWY, m_dofWZ;
    int        m_dofWXP, m_dofWYP, m_dofWZP;
    int        m_dofEF;
    
    DECLARE_FECORE_CLASS();
};

#endif /* FEFluidPResistanceBC_hpp */
