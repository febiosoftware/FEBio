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
#include <FECore/FENodeSetConstraint.h>
#include <FECore/FENodeSet.h>
#include "febiomech_api.h"

class FEBIOMECH_API FEAzimuthConstraint : public FENodeSetConstraint
{
public:
	FEAzimuthConstraint(FEModel* fem);

	//! initialization
	bool Init() override;

	//! Calculate the constraint force
	void LoadVector(FEGlobalVector& R, const FETimeInfo& tp) override;

	//! calculate the constraint stiffness
	void StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp) override;

	//! augmentations
	bool Augment(int naug, const FETimeInfo& tp) override;

	//! build connectivity for matrix profile
	void BuildMatrixProfile(FEGlobalMatrix& M) override;

	//! return number of equations to be allocated for Lagrange multipliers
	int InitEquations(int neq) override;

	//! Update Lagrange multipliers
	void Update(const std::vector<double>& ui) override;

	//! serialization
	void Serialize(DumpStream& ar) override;

	//! return node set
	FENodeSet* GetNodeSet() override;

private:
	int		m_laugon;		//!< Augmentation flag
	double	m_tol;			//!< augmentation tolerance
	double	m_eps;			//!< penalty factor
	int		m_minaug;		//!< minimum nr of augmentations
	int		m_maxaug;		//!< max nr of augmentations

private:
	vector<vec3d>	m_Lm;	//!< Lagrange multipliers

private:
	vector<int>		m_eq;		//!< equation nr of LM
	FEDofList	m_dofU;			//!< dof list
	FENodeSet	m_nodeSet;		//!< the node set to apply this constraint to

	DECLARE_FECORE_CLASS();
};
