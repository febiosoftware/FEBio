//
//  FEFluidResistanceBC.hpp
//  FEBioFluid
//
//  Created by Gerard Ateshian on 9/28/16.
//  Copyright © 2016 febio.org. All rights reserved.
//

#ifndef FEFluidResistanceBC_hpp
#define FEFluidResistanceBC_hpp

#include "FECore/FESurfaceLoad.h"
#include <FECore/FESurfaceMap.h>

//-----------------------------------------------------------------------------
//! FEFluidResistanceBC is a fluid surface that has a normal
//! pressure proportional to the flow rate (resistance).
//!
class FEFluidResistanceBC : public FESurfaceLoad
{
public:
    //! constructor
    FEFluidResistanceBC(FEModel* pfem);
    
    //! Set the surface to apply the load to
    void SetSurface(FESurface* ps) override;
    
    //! calculate traction stiffness (there is none)
    void StiffnessMatrix(const FETimeInfo& tp, FESolver* psolver) override {}
    
    //! calculate residual
    void Residual(const FETimeInfo& tp, FEGlobalVector& R) override { m_alpha = tp.alpha; m_alphaf = tp.alphaf; }
    
    //! mark the dilatation
    void MarkDilatation();
    
    //! set the dilatation
    void SetDilatation();
    
    //! evaluate flow rate
    double FlowRate();
    
    //! initialize
    bool Init() override;

private:
    double			m_R;	//!< flow resistance
    double          m_k;    //!< fluid bulk modulus
    double          m_alpha;
    double          m_alphaf;
    double          m_p0;   //!< fluid pressure offset
    
    int		m_dofWX, m_dofWY, m_dofWZ;
    int		m_dofWXP, m_dofWYP, m_dofWZP;
    int		m_dofEF;
    
    DECLARE_PARAMETER_LIST();
};

#endif /* FEFluidResistanceBC_hpp */
