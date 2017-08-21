#pragma once
#include <FECore/FENLConstraint.h>

//-----------------------------------------------------------------------------
class FERigidBody;

//-----------------------------------------------------------------------------
//! This is a virtual class for all rigid connectors, including
//! spherical, revolute, prismatic and cylindrical joints, as well
//! as springs and dampers that connect rigid bodies.

class FERigidConnector : public FENLConstraint
{
public:
    //! constructor
    FERigidConnector(FEModel* pfem);
    
    //! destructor
    virtual ~FERigidConnector();

	//! initialization
	bool Init();
    
    int GetConnectorID() { return m_nID; }

	//! build connectivity for matrix profile
	void BuildMatrixProfile(FEGlobalMatrix& M);

	//! serialization
	void Serialize(DumpStream& ar);
    
public:
    int	m_nRBa;		//!< rigid body A that the connector connects
    int	m_nRBb;		//!< rigid body B that the connector connects
    
    vec3d	m_F;	//! constraining force
    vec3d	m_M;	//! constraining moment
    
protected:
    int		m_nID;		//!< ID of rigid connector
	bool	m_binit;	//!< initialization flag

	FERigidBody*	m_rbA;
	FERigidBody*	m_rbB;
    
    static int	m_ncount;	//!< used to create unique ID's for the nonlinear constraints
    
    DECLARE_PARAMETER_LIST();
};
