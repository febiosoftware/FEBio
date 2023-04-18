/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

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
#include <FECore/FEModel.h>
#include "febiomech_api.h"

//---------------------------------------------------------------------------------------
class FERigidSystem;
class FERigidBody;
class FERigidPrescribedBC;
class FERigidFixedBC;
class FERigidIC;
class FERigidNodeSet;

//---------------------------------------------------------------------------------------
// This class extends the basic FEModel class by adding a rigid body system
class FEBIOMECH_API FEMechModel : public FEModel
{
public:
	FEMechModel();

	// clear all model data
	void Clear() override;

	// model activation
	void Activate() override;

	// TODO: temporary construction. Would like to call Activate
	void Reactivate() override;

	// reset
	bool Reset() override;

	//! Initialize shells
	void InitShells() override;

	// find a parameter value
	FEParamValue GetParameterValue(const ParamString& param) override;

	//! serialize data for restarts
	void SerializeGeometry(DumpStream& ar) override;

	//! Build the matrix profile for this model
	void BuildMatrixProfile(FEGlobalMatrix& G, bool breset) override;

	// update rigid part of mesh
	void UpdateRigidMesh();

public:
	// get the rigid system
	FERigidSystem* GetRigidSystem();

	// initialize the rigid system
	bool InitRigidSystem() override;

	// number of rigid bodies
	int RigidBodies() const;

	// get a rigid body
	FERigidBody* GetRigidBody(int n);

	// find a rigid body from a material ID
	int FindRigidbodyFromMaterialID(int matId);

	// return number or rigid prescribed BCs
	int RigidPrescribedBCs() const;

	// return the rigid prescribed displacement
	FERigidPrescribedBC* GetRigidPrescribedBC(int i);

	// add a rigid presribed BC
	void AddRigidPrescribedBC(FERigidPrescribedBC* pDC);

	// add a rigid fixed BC
	void AddRigidFixedBC(FERigidFixedBC* pBC);

	// add a rigid initial condition
	void AddRigidInitialCondition(FERigidIC* pIC);

private:
	FERigidSystem*	m_prs;

	DECLARE_FECORE_CLASS();
};
