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
#include <FECore/FESolidElement.h>
#include <FECore/FEGlobalVector.h>
#include <FECore/FEDofList.h>

//-----------------------------------------------------------------------------
//! This class implements a point constraint. That is, it forces a node of a 
//! mesh in the same relative position with respect to the element in which the 
//! node is located.
class FEPointConstraint : public FENLConstraint
{
public:
	//! constructor
	FEPointConstraint(FEModel* pfem);

	//! initialize data
	bool Init() override;

	//! serialize
	void Serialize(DumpStream& ar) override;

	//! Calculate the constraint force
	void LoadVector(FEGlobalVector& R, const FETimeInfo& tp) override;

	//! calculate the constraint stiffness
	void StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp) override;

	//! augmentations \todo implement this
	bool Augment(int naug, const FETimeInfo& tp) override { return true; }

	//! build connectivity for matrix profile
	void BuildMatrixProfile(FEGlobalMatrix& M) override;

public:
	double		m_eps;		//!< penalty parameter
	int			m_node_id;	//!< id of node

public:
	int					m_node;		//!< node to which the constraint is applied
	FESolidElement*		m_pel;		//!< element in which the node is located.
	double				m_rs[3];	//!< natural coordinates in element m_pel

	FEDofList	m_dofU;

	// declare parameter list
	DECLARE_FECORE_CLASS();
};
