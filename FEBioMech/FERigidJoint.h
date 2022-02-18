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
#include "FECore/vec3d.h"
#include "FERigidConnector.h"

//-----------------------------------------------------------------------------
//! The FERigidJoint class implements a rigid joint. The rigid joint allows the
//! user to connect two rigid bodies at a point in space.

class FERigidJoint : public FERigidConnector
{
public:
	//! constructor
	FERigidJoint(FEModel* pfem);

	//! destructor
	~FERigidJoint();

	//! initialization
	bool Init() override;

	// allocate equations
	int InitEquations(int neq) override;

	//! calculates the joint forces
	void LoadVector(FEGlobalVector& R, const FETimeInfo& tp) override;

	//! calculates the joint stiffness
	void StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp) override;

	//! calculate Lagrangian augmentation
	bool Augment(int naug, const FETimeInfo& tp) override;

	//! serialize data to archive
	void Serialize(DumpStream& ar) override;

	//! update state
	void Update() override;

	//! Reset data
	void Reset() override;

protected:
	void UnpackLM(vector<int>& lm);

	// Build the matrix profile
	void BuildMatrixProfile(FEGlobalMatrix& M) override;

	void Update(const std::vector<double>& ui) override;

public:
	vec3d	m_q0;		//! initial position of joint
	vec3d	m_qa0;	//! initial relative position vector of joint w.r.t. A
	vec3d	m_qb0;	//! initial relative position vector of joint w.r.t. B

	vec3d	m_F;		//!< constraining force
	vec3d	m_L;		//!< Lagrange multiplier
	double	m_eps;		//!< penalty factor
	double	m_atol;		//!< augmented Lagrangian tolerance
	int		m_laugon;	//!< enforcement method

protected:
	int		m_nID;	//!< ID of rigid joint

	vector<int>		m_LM;	// Lagrange multiplier equation numbers

	DECLARE_FECORE_CLASS();
};
