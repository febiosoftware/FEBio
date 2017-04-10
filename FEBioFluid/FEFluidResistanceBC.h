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
    void SetSurface(FESurface* ps);
    
    //! calculate traction stiffness (there is none)
    void StiffnessMatrix(const FETimeInfo& tp, FESolver* psolver) {}
    
    //! calculate residual
    void Residual(const FETimeInfo& tp, FEGlobalVector& R) { m_alpha = tp.alpha; }
    
    //! mark the dilatation
    void MarkDilatation();
    
    //! set the dilatation
    void SetDilatation();
    
    //! evaluate flow rate
    double FlowRate();
    
    //! initialize
    bool Init();

private:
    double			m_R;	//!< flow resistance
    double          m_k;    //!< fluid bulk modulus
    double          m_alpha;
    
    int		m_dofVX;
    int		m_dofVY;
    int		m_dofVZ;
    int		m_dofE;
    
    DECLARE_PARAMETER_LIST();
};

#endif /* FEFluidResistanceBC_hpp */
