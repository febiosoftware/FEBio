//
//  FEFluidNormalVelocity.hpp
//  FEBioFluid
//
//  Created by Gerard Ateshian on 3/2/17.
//  Copyright © 2017 febio.org. All rights reserved.
//

#ifndef FEFluidNormalVelocity_hpp
#define FEFluidNormalVelocity_hpp

#include "FECore/FESurfaceLoad.h"
#include <FECore/FESurfaceMap.h>

//-----------------------------------------------------------------------------
//! FEFluidNormalVelocity is a fluid surface that has a normal
//! velocity prescribed on it.
//!
class FEFluidNormalVelocity : public FESurfaceLoad
{
public:
    //! constructor
    FEFluidNormalVelocity(FEModel* pfem);
    
    //! Set the surface to apply the load to
    void SetSurface(FESurface* ps);
    
    //! calculate traction stiffness (there is none)
    void StiffnessMatrix(const FETimeInfo& tp, FESolver* psolver) {}
    
    //! calculate residual
    void Residual(const FETimeInfo& tp, FEGlobalVector& R);
    
    //! Unpack surface element data
    void UnpackLM(FEElement& el, vector<int>& lm);
    
private:
    double			m_velocity;	//!< magnitude of traction load
    FESurfaceMap	m_VC;		//!< traction boundary cards
    
    int		m_dofE;
    
    DECLARE_PARAMETER_LIST();
};

#endif /* FEFluidNormalVelocity_hpp */
