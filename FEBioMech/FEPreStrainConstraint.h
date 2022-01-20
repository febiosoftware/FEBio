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
#include <FECore/FESolidDomain.h>
#include "FEPreStrainElastic.h"
#include "FEPreStrainUncoupledElastic.h"

//-----------------------------------------------------------------------------
class FEPreStrainConstraint : public FENLConstraint
{
public:
	FEPreStrainConstraint(FEModel* pfem);

	// This function must be overloaded by derived classes
	virtual mat3d UpdateFc(const mat3d& F, const mat3d& Fc_prev, FEMaterialPoint& mp, FEPrestrainMaterial* pmat) = 0;

public:
	bool Init() override;
	bool Augment(int naug, const FETimeInfo& tp) override;

	void LoadVector(FEGlobalVector& R, const FETimeInfo& tp) override {};
	void StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp) override {};
	virtual void Update(const FETimeInfo& tp) {}

	void BuildMatrixProfile(FEGlobalMatrix& M) override {}

private:
	bool Augment(FESolidDomain* pdom, int n, int naug);

public:
	bool	m_laugon;	//!< augmented Lagrangian flag
	double	m_tol;		//!< convergence tolerance
	int		m_naugmin;	//!< minimum number of iterations
	int		m_naugmax;	//!< maximum number of iterations

	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
// Implements a constraint that essentially eliminates the distortion
class FEGPAConstraint : public FEPreStrainConstraint
{
public:
	FEGPAConstraint(FEModel* pfem) : FEPreStrainConstraint(pfem){}

	mat3d UpdateFc(const mat3d& F, const mat3d& Fc_prev, FEMaterialPoint& mp, FEPrestrainMaterial* pmat) override;
};

//-----------------------------------------------------------------------------
// enforces just the fiber stretch
class FEInSituStretchConstraint: public FEPreStrainConstraint
{
public:
	FEInSituStretchConstraint(FEModel* pfem);
	mat3d UpdateFc(const mat3d& F, const mat3d& Fc_prev, FEMaterialPoint& mp, FEPrestrainMaterial* pmat) override;

private:
	double	m_max_stretch;
	bool	m_biso;

	DECLARE_FECORE_CLASS();
};
