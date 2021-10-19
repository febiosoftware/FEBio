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
#include <FECore/vec3d.h>
#include "FERigidConnector.h"
#include "febiomech_api.h"

//-----------------------------------------------------------------------------
class FEBIOMECH_API FEGenericRigidJoint : public FERigidConnector
{
public:
	FEGenericRigidJoint(FEModel* fem);

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

	//! initial position 
	vec3d InitialPosition() const;

	//! current position
	vec3d Position() const;

protected:
	void UnpackLM(vector<int>& lm);

	// Build the matrix profile
	void BuildMatrixProfile(FEGlobalMatrix& M) override;

	void Update(const std::vector<double>& Ui, const std::vector<double>& ui) override;
	void UpdateIncrements(std::vector<double>& Ui, const std::vector<double>& ui) override;

	void PrepStep() override;

protected:
	//! stiffness matrix for penalty and augmented Lagrangian
	void StiffnessMatrixAL(FELinearSystem& LS, const FETimeInfo& tp);

	//! stiffness matrix for Lagrange Multipliers
	void StiffnessMatrixLM(FELinearSystem& LS, const FETimeInfo& tp);

private:
	int		m_laugon;
	double	m_eps;
	double	m_tol;
	vec3d	m_q0;
	bool	m_bsymm;

private:
	vec3d	m_qa0;
	vec3d	m_qb0;

	bool	m_bc[6];
	vec3d	m_v0[3];
	double	m_dc[6];

	int		m_EQ[6];
	double	m_Lm[6], m_Lmp[6];

	DECLARE_FECORE_CLASS();
};
