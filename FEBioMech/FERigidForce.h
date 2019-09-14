/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in 
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

//-----------------------------------------------------------------------------
//! an axial force between two rigid bodies
class FERigidAxialForce : public FEModelLoad
{
public:
	//! constructor
	FERigidAxialForce(FEModel* pfem);

	//! initialization
	bool Init() override;

	//! Serialization
	void Serialize(DumpStream& ar) override;

	//! Residual
	void LoadVector(FEGlobalVector& R, const FETimeInfo& tp) override;

	//! Stiffness matrix
	void StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp) override;

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
class FERigidBodyForce : public FEModelLoad
{
public:
	enum { RAMP, TARGET };	// values for m_ntype

public:
	FERigidBodyForce(FEModel* pfem);

	//! Activation
	void Activate() override;

	//! initialization
	bool Init() override;

	//! get the current force value
	double Value();

	//! Serialization
	void Serialize(DumpStream& ar) override;

	//! forces
	void LoadVector(FEGlobalVector& R, const FETimeInfo& tp) override;

	//! Stiffness matrix
	void StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp) override;

public:
	void SetType(int ntype) { m_ntype = ntype; }

	void SetID(int nid) { m_id = nid; }

	void SetBC(int bc) { m_bc = bc; }

	void SetFollowFlag(bool b) { m_bfollow = b; }

	void SetForce(double f) { m_force = f; }

private:
	int		m_ntype;		//!< type of force (0=loadcurve, 1=target)
	int		m_id;			//!< rigid body id
	int		m_bc;			//!< force direction
	double	m_force;		//!< applied force
	double	m_trg;			//!< target force for target case
	bool	m_bfollow;		//!< follower force if true

	DECLARE_FECORE_CLASS();
};
