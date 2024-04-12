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
#include <FECore/FESurfaceConstraint.h>
#include <FECore/FEEdge.h>
#include "FEContactSurface.h"
#include <set>

class FEE2SCPPoint : public FELineMaterialPoint
{
public:
	double	m_gap = 0.0;
	double	m_Ln = 0.0;
	vec3d	m_tc;

	void Serialize(DumpStream& ar) override
	{
		FELineMaterialPoint::Serialize(ar);
		ar& m_tc & m_gap;
	}
};

class FEE2SCPSurface : public FESurface
{
public:
	FEE2SCPSurface(FEModel* fem);

	void Update();
};

class FEE2SCPEdge : public FEEdge
{
public:
	FEE2SCPEdge(FEModel* fem);

	FEMaterialPoint* CreateMaterialPoint() override;

	void Update();

	bool Create(FESegmentSet& eset) override;
};

class FEEdgeToSurfaceContactPotential : public FESurfaceConstraint
{
public:
	FEEdgeToSurfaceContactPotential(FEModel* fem);

public:
	//! return the surface
	FESurface* GetSurface() override;

	// Build the matrix profile
	void BuildMatrixProfile(FEGlobalMatrix& M) override;

	// update
	void Update() override;

	// init
	bool Init() override;

	// serialization
	void Serialize(DumpStream& ar) override;

public:
	// The LoadVector function evaluates the "forces" that contribute to the residual of the system
	void LoadVector(FEGlobalVector& R, const FETimeInfo& tp) override;

	// Evaluates the contriubtion to the stiffness matrix
	void StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp) override;

protected:
	void ElementForce(FELineElement& el1, FESurfaceElement& el2, std::vector<double>& fe);
	void ElementStiffness(FELineElement& el1, FESurfaceElement& el2, matrix& ke);

	double PotentialDerive(double r);
	double PotentialDerive2(double r);

protected:
	FEE2SCPEdge	m_edge;
	FEE2SCPSurface	m_surf;

protected:
	double	m_kc;
	double	m_p;
	double	m_Rin;
	double	m_Rout;
	double	m_Rmin;
	double	m_wtol;

	double	m_c1, m_c2;

	std::vector< std::set<FESurfaceElement*> >	m_activeElements;

	DECLARE_FECORE_CLASS();
};
