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

	//! Clear
	void Clear();

	//! Initialization
	bool Init();

	//! Reset data
	bool Reset();

	//! place data on stream for running restarts
	void ShallowCopy(DumpStream& dmp, bool bsave);

	// Find a parameter (from a rigid material index)
	double* FindParameter(int nmat, ParamString& sz, int index);

protected:
	bool CreateObjects();

private:
	FEModel&					m_fem;	//!< the FE model this system is attached to
	std::vector<FERigidBody*>	m_RB;	//!< the list of rigid bodies in this system
};
