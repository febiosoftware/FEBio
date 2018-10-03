#pragma once
#include "FERigidConnector.h"

//-----------------------------------------------------------------------------
//! The FERigidAngularDamper class implements an angular damper that connects
//! two rigid bodies.

class FERigidAngularDamper : public FERigidConnector
{
public:
    //! constructor
    FERigidAngularDamper(FEModel* pfem);
    
    //! destructor
    ~FERigidAngularDamper() {}
    
    //! initialization
    bool Init() override;
    
    //! calculates the joint forces
    void Residual(FEGlobalVector& R, const FETimeInfo& tp) override;
    
    //! calculates the joint stiffness
    void StiffnessMatrix(FESolver* psolver, const FETimeInfo& tp) override;
    
    //! calculate Lagrangian augmentation
    bool Augment(int naug, const FETimeInfo& tp) override;
    
    //! serialize data to archive
    void Serialize(DumpStream& ar) override;
    
    //! update state
    void Update(int niter, const FETimeInfo& tp) override;
    
    //! Reset data
    void Reset() override;
    
public: // parameters
    double	m_c;        //! damping constant
    
    DECLARE_FECORE_CLASS();
};
