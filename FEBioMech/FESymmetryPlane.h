#pragma once
#include "FEAugLagLinearConstraint.h"
#include <FECore/FESurface.h>

//-----------------------------------------------------------------------------
//! The FESymmetryPlane class implements a symmetry plane
//! as a linear constraint on the components of the solid displacement.

class FESymmetryPlane : public FELinearConstraintSet
{
public:
    //! constructor
    FESymmetryPlane(FEModel* pfem);
    
    //! destructor
    ~FESymmetryPlane() {}
    
    //! Activation
    void Activate() override;
    
    //! initialization
    bool Init() override;
    
    //! Get the surface
    FESurface* GetSurface() override { return &m_surf; }
    
protected:
    FESurface    m_surf;
};
