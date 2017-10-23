#pragma once
#include "FECore/vec3d.h"
#include "FERigidConnector.h"

//-----------------------------------------------------------------------------
//! The FERigidDamper class implements a linear damper that connects
//! two rigid bodies at arbitrary points (not necessarily nodes).
//! TODO: This inherits from FENLConstraint, which is not the appropriate base class
class FERigidDamper : public FERigidConnector
{
public:
    //! constructor
    FERigidDamper(FEModel* pfem);
    
    //! destructor
    ~FERigidDamper() {}
    
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
	void Update(int niter, const FETimeInfo& tp);
    
    //! Reset data
    void Reset();
    
public: // parameters
    double	m_c;        //! damping constant
    vec3d	m_a0;       //! initial absolute position vector of spring on body A
    vec3d	m_b0;       //! initial absolute position vector of spring on body B

protected:
    vec3d	m_qa0;      //! initial relative position vector of spring on body A
    vec3d	m_qb0;      //! initial relative position vector of spring on body B
    
    DECLARE_PARAMETER_LIST();
};
