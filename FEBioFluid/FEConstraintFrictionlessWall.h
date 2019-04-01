#pragma once
#include <FEBioMech/FEAugLagLinearConstraint.h>
#include <FECore/FESurface.h>
#include "febiofluid_api.h"

//-----------------------------------------------------------------------------
//! The FEConstraintFrictionlessWall class implements a frictionless fluid wall
//! as a linear constraint on the components of the fluid velocity.

class FEBIOFLUID_API FEConstraintFrictionlessWall : public FELinearConstraintSet
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
    FESurface* GetSurface() { return &m_surf; }
    
protected:
    FESurface	m_surf;
};
