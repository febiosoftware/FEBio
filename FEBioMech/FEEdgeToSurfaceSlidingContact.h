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

//=======================================================================================
class FEE2SSlidingContactPoint : public FELineMaterialPoint
{
public:
	double	m_gap = 0.0;
	double	m_Ln = 0.0;	  //!< Lagrange multipliers for contact pressure
	double	m_Lm = 0.0;	  //!< Lagrange multipliers for contact pressure
	vec3d	m_nu;	  //!< secondary surface normal at primary surface node
	vec2d	m_rs;	  //!< natural coordinates of primary surface projection on secondary surface element
	FESurfaceElement* m_pme = nullptr;

	void Serialize(DumpStream& ar) override
	{
		FELineMaterialPoint::Serialize(ar);
		ar & m_gap & m_Ln & m_Lm & m_nu & m_rs;
	}

	void Init() override
	{
		FELineMaterialPoint::Init();
		m_gap = 0.0;
		m_Ln = 0.0;
		m_Lm = 0.0;
		m_nu = vec3d(0, 0, 0);
		m_rs = vec2d(0, 0);
		m_pme = nullptr;
	}
};

//=======================================================================================
class FEEdgeToSurfaceSlidingContactSurface : public FEContactSurface
{
public:
	FEEdgeToSurfaceSlidingContactSurface(FEModel* fem);

	void Update();
};

//=======================================================================================
class FEEdgeToSurfaceSlidingContactEdge : public FEEdge
{
public:
	FEEdgeToSurfaceSlidingContactEdge(FEModel* fem);

	FEMaterialPoint* CreateMaterialPoint() override;

	bool Init() override;

	void Update();

	bool Create(FESegmentSet& eset) override;

	void UnpackLM(FELineElement& el, vector<int>& lm);

	void Serialize(DumpStream& ar);

public:
	int m_dofX = -1;
	int m_dofY = -1;
	int m_dofZ = -1;

	vector<FEE2SSlidingContactPoint> m_points;	//!< sliding contact surface data
};

//=======================================================================================
class FEEdgeToSurfaceSlidingContact : public FESurfaceConstraint
{
public:
	FEEdgeToSurfaceSlidingContact(FEModel* fem);

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

	void ProjectSurface(bool bupseg, bool bmove = false);

	void ContactNodalForce(int m, FESurfaceElement& mel, vector<double>& fe);
	void ContactNodalStiffness(int m, FESurfaceElement& mel, matrix& ke);

protected:
	FEEdgeToSurfaceSlidingContactEdge	m_edge;
	FEEdgeToSurfaceSlidingContactSurface	m_surf;

protected:
	double	m_atol;
	double	m_eps;
	int		m_naugmin;
	int		m_naugmax;
	double	m_stol;
	int		m_nsegup;
	double	m_sradius;

	bool m_bfirst;

	DECLARE_FECORE_CLASS();
};
