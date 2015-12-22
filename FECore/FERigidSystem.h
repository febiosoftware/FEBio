#pragma once
#include "FECoreBase.h"
#include "FERigidBody.h"
#include <vector>

//-----------------------------------------------------------------------------
// forward declarations
class FEModel;

//-----------------------------------------------------------------------------
//! The FERigidSystem class manages all rigid body paraphernalia.
class FERigidSystem
{
public:
	//! constructor
	FERigidSystem(FEModel* pfem);

	//! Add a rigid body
	void AddRigidBody(FERigidBody* prb);

	//! return the number of rigid bodies
	int Objects() const;

	//! get a rigid body
	FERigidBody* Object(int i);

	//! Activate
	void Activate();

	//! Serialize
	void Serialize(DumpFile& ar);

	//! Clear
	void Clear();

	//! Initialization
	bool Init();

	//! Reset data
	bool Reset();

	// find a model component from its ID
	FEModelComponent* FindModelComponent(int nid);

	//! place data on stream for running restarts
	void ShallowCopy(DumpStream& dmp, bool bsave);

	// Find a parameter (from a rigid material index)
	double* FindParameter(int nmat, ParamString& sz, int index);

	// update the mesh geometry
	void UpdateMesh();

public:
	// rigid nodes
	int RigidNodes() { return (int) m_RN.size(); }
	FERigidNode* RigidNode(int i) { return m_RN[i]; }
	void AddRigidNode(FERigidNode* prn) { m_RN.push_back(prn); }

protected:
	bool CreateObjects();

public:
	// Boundary/Initial conditions for rigid bodies
	// TODO: I'd like to do something different with this. Perhaps place them in the BC or in some constraint section.
	vector<FERigidNode*>				m_RN;	//!< rigid nodes
	vector<FERigidBodyFixedBC*>			m_RBC;	//!< rigid body fixed
	vector<FERigidBodyDisplacement*>	m_RDC;	//!< rigid body displacements
	vector<FERigidBodyVelocity*>		m_RBV;	//!< rigid body initial velocities
	vector<FERigidBodyAngularVelocity*>	m_RBW;	//!< rigid body initial angular velocities

private:
	FEModel&					m_fem;	//!< the FE model this system is attached to
	std::vector<FERigidBody*>	m_RB;	//!< the list of rigid bodies in this system
};
