/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, Columbia University, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/

#pragma once
#include <FECore/FECoreBase.h>
#include <vector>
#include "febiomech_api.h"

//-----------------------------------------------------------------------------
// forward declarations
class FEModel;
class FERigidBody;
class FEModelComponent;
class FERigidNodeSet;
class FERigidBodyFixedBC;
class FERigidBodyDisplacement;
class FERigidBodyVelocity;
class FERigidBodyAngularVelocity;
class FERigidSurface;
class FEGlobalMatrix;

//-----------------------------------------------------------------------------
//! The FERigidSystem class manages all rigid body paraphernalia.
class FEBIOMECH_API FERigidSystem
{
public:
	//! constructor
	FERigidSystem(FEModel* pfem);

	//! destructor
	~FERigidSystem();

	//! Add a rigid body
	void AddRigidBody(FERigidBody* prb);

	//! return the number of rigid bodies
	int Objects() const;

	//! get a rigid body
	FERigidBody* Object(int i);

	//! Activate
	void Activate();

	//! Clear
	void Clear();

	//! Initialization
	bool Init();

	//! Reset data
	bool Reset();

	//! place data on stream for restarts
	void Serialize(DumpStream& dmp);

	// Find a parameter (from a rigid material index)
	double* FindParameter(int nmat, ParamString& sz, int index);

	FEParamValue GetParameterValue(const ParamString& paramString);

	// update the mesh geometry
	void UpdateMesh();

	// build the matrix profile for the rigid system
	void BuildMatrixProfile(FEGlobalMatrix& G);

	// find a rigid body from a material ID
	int FindRigidbodyFromMaterialID(int matId);

public:
	int RigidNodeSets() { return (int) m_RN.size(); }
	FERigidNodeSet* RigidNodeSet(int i) { return m_RN[i]; }
	void AddRigidNodeSet(FERigidNodeSet* prn) { m_RN.push_back(prn); }

	int FixedBCs() { return (int) m_RBC.size(); }
	FERigidBodyFixedBC* FixedBC(int i) { return m_RBC[i]; }
	void AddFixedBC(FERigidBodyFixedBC* pbc) { m_RBC.push_back(pbc); }

	int PrescribedBCs() { return (int) m_RDC.size(); }
	FERigidBodyDisplacement* PrescribedBC(int i) { return m_RDC[i]; }
	void AddPrescribedBC(FERigidBodyDisplacement* pdc) { m_RDC.push_back(pdc); }

	int InitialVelocities() { return (int) m_RBV.size(); }
	FERigidBodyVelocity* InitialVelocity(int i) { return m_RBV[i]; }
	void AddInitialVelocity(FERigidBodyVelocity* pbv) { m_RBV.push_back(pbv); }

	int InitialAngularVelocities() { return (int) m_RBW.size(); }
	FERigidBodyAngularVelocity* InitialAngularVelocity(int i) { return m_RBW[i]; }
	void AddInitialAngularVelocity(FERigidBodyAngularVelocity* pbw) { m_RBW.push_back(pbw); }

	std::vector<FERigidBody*>& RigidBodyList();

public:
	void AddRigidSurface(FERigidSurface* prs);

	FERigidSurface* FindRigidSurface(const std::string& name);

protected:
	bool CreateObjects();

protected:
	// Boundary/Initial conditions for rigid bodies
	// TODO: I'd like to do something different with this. Perhaps place them in the BC or in some constraint section.
	vector<FERigidNodeSet*>				m_RN;	//!< rigid node sets
	vector<FERigidBodyFixedBC*>			m_RBC;	//!< rigid body fixed
	vector<FERigidBodyDisplacement*>	m_RDC;	//!< rigid body displacements
	vector<FERigidBodyVelocity*>		m_RBV;	//!< rigid body initial velocities
	vector<FERigidBodyAngularVelocity*>	m_RBW;	//!< rigid body initial angular velocities

private:
	FEModel&					m_fem;	//!< the FE model this system is attached to
	vector<FERigidBody*>	m_RB;	//!< the list of rigid bodies in this system
	vector<FERigidSurface*>	m_RS;	//!< rigid surfaces
};
