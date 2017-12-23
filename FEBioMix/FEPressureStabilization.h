//
//  FEPressureStabilization.hpp
//  FEBioMix
//
//  Created by Gerard Ateshian on 8/8/17.
//  Copyright © 2017 febio.org. All rights reserved.
//

#ifndef FEPressureStabilization_hpp
#define FEPressureStabilization_hpp

#include "FECore/FESurfaceLoad.h"
#include <FECore/FESurfaceMap.h>

//-----------------------------------------------------------------------------
//! This pseudo-surface load is used to calculate the pressure stabilization
//! time constant based on the properties of elements under that surface
//!
class FEPressureStabilization : public FESurfaceLoad
{
public:
    //! constructor
    FEPressureStabilization(FEModel* pfem);
    
    //! Set the surface to apply the load to
    void SetSurface(FESurface* ps) override;
    
    //! calculate pressure stiffness
    void StiffnessMatrix(const FETimeInfo& tp, FESolver* psolver) override {}
    
    //! calculate residual
    void Residual(const FETimeInfo& tp, FEGlobalVector& R) override {}
    
    //! initialize
    bool Init() override;
    
protected:
    double TimeConstant(FESurfaceElement& el, FESurface& s);
    
protected:
    bool	m_bstab;		//!< flag for calculating stabilization constant
    
    DECLARE_PARAMETER_LIST();
};

#endif /* FEPressureStabilization_hpp */
