#pragma once
#include "FECore/vec3d.h"
#include "FERigidConnector.h"

//-----------------------------------------------------------------------------
//! The FERigidContractileForce class implements a contractile force between
//! arbitrary points (not necessarily nodes) on two rigid bodies.

class FERigidContractileForce : public FERigidConnector
{
public:
    //! constructor
    FERigidContractileForce(FEModel* pfem);
    
    //! destructor
    ~FERigidContractileForce() {}
    
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
    double	m_f0;       //! contractile force
    vec3d	m_a0;       //! initial absolute position vector of insertion on body A
    vec3d	m_b0;       //! initial absolute position vector of insertion on body B

protected:
    vec3d	m_qa0;      //! initial relative position vector of insertion on body A
    vec3d	m_qb0;      //! initial relative position vector of insertion on body B
    
    DECLARE_FECORE_CLASS();
};
