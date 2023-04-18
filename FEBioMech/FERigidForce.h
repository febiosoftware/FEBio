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
#include "FECore/FEModelLoad.h"
#include "febiomech_api.h"

//-----------------------------------------------------------------------------
class FEBIOMECH_API FERigidLoad : public FEModelLoad
{
public:
	FERigidLoad(FEModel* fem) : FEModelLoad(fem) {}

	FECORE_BASE_CLASS(FERigidLoad)
};

//-----------------------------------------------------------------------------
//! an axial force between two rigid bodies
class FEBIOMECH_API FERigidAxialForce : public FERigidLoad
{
public:
	//! constructor
	FERigidAxialForce(FEModel* pfem);

	//! initialization
	bool Init() override;

	//! Serialization
	void Serialize(DumpStream& ar) override;

	//! Residual
	void LoadVector(FEGlobalVector& R) override;

	//! Stiffness matrix
	void StiffnessMatrix(FELinearSystem& LS) override;

public:
	int		m_ida, m_idb;		//!< rigid body ID's
	vec3d	m_ra0, m_rb0;		//!< coordinates of attachements in reference state
	double	m_s;				//!< scale factor
	bool	m_brelative;		//!< if active, the ra0 and rb0 are relative w.r.t. the COM

	DECLARE_FECORE_CLASS();
};


//-----------------------------------------------------------------------------
//! rigid body force
//! TODO: I'd like to split this class into two classes: one handling the case
//!       were the force is const, and one where the force is a follower force.
//!       Perhaps I can derive the const force from FENodalLoad since it applies
//!       a force directly to the rigid "node".
class FEBIOMECH_API FERigidBodyForce : public FERigidLoad
{
public:
	enum { FORCE_LOAD, FORCE_FOLLOW, FORCE_TARGET };	// values for m_ntype

public:
	FERigidBodyForce(FEModel* pfem);

	//! Activation
	void Activate() override;

	//! initialization
	bool Init() override;

	//! Serialization
	void Serialize(DumpStream& ar) override;

	//! forces
	void LoadVector(FEGlobalVector& R) override;

	//! Stiffness matrix
	void StiffnessMatrix(FELinearSystem& LS) override;

public:
	void SetRigidMaterialID(int nid);

	void SetDOF(int bc);

	void SetLoadType(int loadType);

	void SetForce(double f);

private:
	int		m_rigidMat;		//!< rigid body material id
	int		m_dof;			//!< force direction
	bool	m_brelative;	//!< relative flag

	int		m_ntype;		//!< type of force
	double	m_force;		//!< applied force
	double	m_force0;		//!< initial force at activation
	int		m_rid;

	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
//! rigid body moment
class FEBIOMECH_API FERigidBodyMoment : public FERigidLoad
{
public:
	enum { MOMENT_LOAD, MOMENT_FOLLOW, MOMENT_TARGET };	// values for m_ntype

public:
	FERigidBodyMoment(FEModel* pfem);

	//! Activation
	void Activate() override;

	//! initialization
	bool Init() override;

	//! Serialization
	void Serialize(DumpStream& ar) override;

	//! forces
	void LoadVector(FEGlobalVector& R) override;

	//! Stiffness matrix
	void StiffnessMatrix(FELinearSystem& LS) override;

public:
	void SetRigidMaterialID(int nid);

	void SetDOF(int bc);

	void SetValue(double f);

private:
	int		m_rigidMat;		//!< rigid body material id
	int		m_dof;			//!< force direction
	bool	m_brelative;	//!< relative flag
	double	m_value;		//!< applied moment
	int		m_ntype;		//!< type of moment

	double	m_value0;		//!< initial moment at activation (used with brelative flag)
	int		m_rid;

	DECLARE_FECORE_CLASS();
};
