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
#include <FECore/FENLConstraint.h>
#include <FECore/FEDofList.h>

//-----------------------------------------------------------------------------
// This class implements a constraint that enforces the distance between two nodes
class FEDistanceConstraint : public FENLConstraint
{
public:
	//! constructor
	FEDistanceConstraint(FEModel* pfem);

	bool Init() override;
	void Activate() override;
	void LoadVector(FEGlobalVector& R, const FETimeInfo& tp) override;
	void StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp) override;
	bool Augment(int naug, const FETimeInfo& tp) override;
	void Serialize(DumpStream& ar) override;

	//! build connectivity for matrix profile
	void BuildMatrixProfile(FEGlobalMatrix& M) override;

	// update state
	void Reset() override;

public:
	double	m_eps;		//!< penalty parameter
	double	m_atol;		//!< augmented Lagrangian tolerance
	bool	m_blaugon;	//!< augmentation flag
	int		m_node[2];	//!< the two nodes that are connected
	int		m_nminaug;	//!< min number of augmentations
	int		m_nmaxaug;	//!< max number of augmentations

	double	m_l0;		//!< reference length
	double	m_Lm;		//!< Lagrange multiplier

	FEDofList	m_dofU;

	DECLARE_FECORE_CLASS();
};
