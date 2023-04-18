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

#include <FECore/vector.h>
#include <FECore/matrix.h>
#include <FECore/FESurfaceConstraint.h>
#include <FECore/FENodeSetConstraint.h>
#include <FECore/FECoreClass.h>
#include "fecore_api.h"
#include <list>

//-----------------------------------------------------------------------------
//! linear constraint enforced using augmented lagrangian

class FECORE_API FEAugLagLinearConstraintDOF : public FECoreClass
{
public:
	FEAugLagLinearConstraintDOF(FEModel* fem);

public:
	int	m_node;		// the node to which this dof belongs to
	int	m_bc;			// the degree of freedom
	double	m_val;	// coefficient value

	DECLARE_FECORE_CLASS();
	FECORE_BASE_CLASS(FEAugLagLinearConstraintDOF);
};

class FECORE_API FEAugLagLinearConstraint : public FECoreClass
{
public:
	typedef std::vector<FEAugLagLinearConstraintDOF*>::iterator Iterator;

public:
	//! constructor
	FEAugLagLinearConstraint(FEModel* fem) : FECoreClass(fem) { m_lam = 0; }

	//! serialize data to archive
	void Serialize(DumpStream& ar);

	void ClearDOFs();

	void AddDOF(int node, int bc, double val);

public:
	std::vector<FEAugLagLinearConstraintDOF*>	m_dof;	//!< list of participating dofs
	double									m_lam;	//!< lagrange multiplier

	DECLARE_FECORE_CLASS();
	FECORE_BASE_CLASS(FEAugLagLinearConstraint);
};

//-----------------------------------------------------------------------------
//! This class manages a group of linear constraints

class FECORE_API FELinearConstraintSet : public FENLConstraint
{
public:
	//! constructor
	FELinearConstraintSet(FEModel* pfem);

	//! add a linear constraint to the list
	void add(FEAugLagLinearConstraint* plc) { m_LC.push_back(plc); }

public:
	//! add the linear constraint contributions to the residual
	void LoadVector(FEGlobalVector& R, const FETimeInfo& tp) override;

	//! add the linear constraint contributions to the stiffness matrix
	void StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp) override;

	//! do the augmentation
	bool Augment(int naug, const FETimeInfo& tp) override;

	//! build connectivity for matrix profile
	void BuildMatrixProfile(FEGlobalMatrix& M) override;

	//! serialization
	void Serialize(DumpStream& ar) override;

protected:
	//! calculate the constraint value
	double constraint(FEAugLagLinearConstraint& LC);

public:
	std::vector<FEAugLagLinearConstraint*>	m_LC;	//!< list of linear constraints

public:
	bool	m_laugon;	//!< augmentation flag
	double	m_tol;	//!< augmentation tolerance
	double	m_eps;	//!< penalty factor
    double  m_rhs;  //!< right-hand-side of linear constraint equation
	int		m_naugmax;	//!< max nr of augmentations
	int		m_naugmin;	//!< min nf of augmentations

	int	m_nID;		//!< ID of manager

	DECLARE_FECORE_CLASS();
};
