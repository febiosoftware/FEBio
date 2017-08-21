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
    bool Init();
    
    //! calculates the joint forces
    void Residual(FEGlobalVector& R, const FETimeInfo& tp);
    
    //! calculates the joint stiffness
    void StiffnessMatrix(FESolver* psolver, const FETimeInfo& tp);
    
    //! calculate Lagrangian augmentation
    bool Augment(int naug, const FETimeInfo& tp);
    
    //! serialize data to archive
    void Serialize(DumpStream& ar);
    
    //! update state
    void Update(const FETimeInfo& tp);
    
    //! Reset data
    void Reset();
    
public: // parameters
    double	m_c;        //! damping constant
    
    DECLARE_PARAMETER_LIST();
};
