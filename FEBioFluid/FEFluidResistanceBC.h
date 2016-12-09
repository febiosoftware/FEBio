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
    void StiffnessMatrix(const FETimeInfo& tp, FESolver* psolver);
    
    //! calculate residual
    void Residual(const FETimeInfo& tp, FEGlobalVector& R);
    
    //! Unpack surface element data
    void UnpackLM(FEElement& el, vector<int>& lm);
    
    //! evaluate flow rate
    double FlowRate();
    
    //! calculate stiffness for an element
    void FluidResistanceStiffness(FESurfaceElement& el, matrix& ke);

private:
    double			m_R;	//!< flow resistance
    
    int		m_dofVX;
    int		m_dofVY;
    int		m_dofVZ;
    int		m_dofE;
    
    DECLARE_PARAMETER_LIST();
};

#endif /* FEFluidResistanceBC_hpp */
