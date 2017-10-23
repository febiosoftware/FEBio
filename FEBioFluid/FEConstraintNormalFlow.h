//
//  FEConstraintNormalFlow.hpp
//  FEBioFluid
//
//  Created by Gerard Ateshian on 9/5/16.
//  Copyright © 2016 febio.org. All rights reserved.
//

#ifndef ConstraintNormalFlow_hpp
#define ConstraintNormalFlow_hpp

#include <FEBioMech/FEAugLagLinearConstraint.h>
#include <FECore/FESurface.h>

//-----------------------------------------------------------------------------
//! The FEConstraintNormalFlow class implements a fluid surface with zero
//! tangential velocity as a linear constraint.

class FEConstraintNormalFlow : public FELinearConstraintSet
{
public:
    //! constructor
    FEConstraintNormalFlow(FEModel* pfem);
    
    //! destructor
    ~FEConstraintNormalFlow() {}
    
    //! Activation
    void Activate();
    
    //! initialization
    bool Init();
    
    //! Get the surface
    FESurface* GetSurface() { return &m_surf; }
    
protected:
    FESurface	m_surf;
};

#endif /* ConstraintNormalFlow_hpp */
