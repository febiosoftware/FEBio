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
#include "FERigidForce.h"
#include <FECore/FECoreClass.h>

class FERigidBody;

class FEBIOMECH_API FERigidCablePoint : public FECoreClass
{
public:
	FERigidCablePoint(FEModel* fem) : FECoreClass(fem) {}

public:
	int		m_rb;	//!< rigid body ID
	vec3d	m_pos;	//!< position of attachment point

	DECLARE_FECORE_CLASS();
	FECORE_BASE_CLASS(FERigidCablePoint);
};

class FERigidCable : public FERigidLoad
{
public:
	FERigidCable(FEModel* fem);

	//! initialization
	bool Init() override;

	//! forces
	void LoadVector(FEGlobalVector& R) override;

	//! Stiffness matrix
	void StiffnessMatrix(FELinearSystem& LS) override;

private:
	void applyRigidForce(FERigidBody& rb, const vec3d& F, const vec3d& d, FEGlobalVector& R);

private:
	double	m_force;		//!< magnitude of force (i.e. tension in cable)
	vec3d	m_forceDir;		//!< direction of force at cable's end
	bool	m_brelative;	//!< positions are defined relative w.r.t. rigid body's COM or not
	std::vector<FERigidCablePoint*>	m_points;

private:
	DECLARE_FECORE_CLASS();
};
