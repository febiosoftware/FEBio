#pragma once
#include <FEBioMech/FEAugLagLinearConstraint.h>
#include <FECore/FESurface.h>
#include "febiofluid_api.h"

//-----------------------------------------------------------------------------
//! The FEConstraintNormalFlow class implements a fluid surface with zero
//! tangential velocity as a linear constraint.

class FEBIOFLUID_API FEConstraintNormalFlow : public FELinearConstraintSet
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
