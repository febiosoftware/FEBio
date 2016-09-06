//
//  FEConstraintFrictionlessWall.hpp
//  FEBioFluid
//
//  Created by Gerard Ateshian on 9/5/16.
//  Copyright © 2016 febio.org. All rights reserved.
//

#ifndef FEConstraintFrictionlessWall_hpp
#define FEConstraintFrictionlessWall_hpp

#include <FEBioMech/FEAugLagLinearConstraint.h>

//-----------------------------------------------------------------------------
//! The FEConstraintFrictionlessWall class implements a frictionless fluid wall
//! as a linear constraint on the components of the fluid velocity.

class FEConstraintFrictionlessWall : public FELinearConstraintSet
{
public:
    //! constructor
    FEConstraintFrictionlessWall(FEModel* pfem);
    
    //! destructor
    ~FEConstraintFrictionlessWall() {}
    
    //! Activation
    void Activate();

    //! initialization
    bool Init();
    
    //! Get the surface
    FESurface* GetSurface(const char* sz) { return &m_surf; }
    
protected:
    FESurface	m_surf;
};

#endif /* FEConstraintFrictionlessWall_hpp */
