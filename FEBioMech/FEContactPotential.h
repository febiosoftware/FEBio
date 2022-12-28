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
#include "FEContactInterface.h"
#include "FEContactSurface.h"
#include <set>

class FEContactPotentialSurface : public FEContactSurface
{
public:
	class Data : public FEContactMaterialPoint
	{
	public:
		vec3d	m_tc;

		void Serialize(DumpStream& ar) override
		{
			FEContactMaterialPoint::Serialize(ar);
			ar & m_tc;
		}
	};

public:
	FEContactPotentialSurface(FEModel* fem);

	void GetContactTraction(int nelem, vec3d& tc) override;

	double GetContactArea() override;

	FEMaterialPoint* CreateMaterialPoint() override;
};

typedef FEContactPotentialSurface::Data FECPContactPoint;

class FEContactPotential : public FEContactInterface
{
public:
	FEContactPotential(FEModel* fem);

	// -- From FESurfacePairConstraint
public:
	//! return the primary surface
	FESurface* GetPrimarySurface() override;

	//! return the secondary surface
	FESurface* GetSecondarySurface() override;

	//! temporary construct to determine if contact interface uses nodal integration rule (or facet)
	bool UseNodalIntegration() override;

	// Build the matrix profile
	void BuildMatrixProfile(FEGlobalMatrix& M) override;

	// update
	void Update() override;

	// init
	bool Init() override;

	// serialization
	void Serialize(DumpStream& ar) override;

	// -- from FEContactInterface
public:
	// The LoadVector function evaluates the "forces" that contribute to the residual of the system
	void LoadVector(FEGlobalVector& R, const FETimeInfo& tp) override;

	// Evaluates the contriubtion to the stiffness matrix
	void StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp) override;

protected:
	void ElementForce(FESurfaceElement& el1, FESurfaceElement& el2, std::vector<double>& fe);
	void ElementStiffness(FESurfaceElement& el1, FESurfaceElement& el2, matrix& ke);

	double PotentialDerive(double r);
	double PotentialDerive2(double r);

	void BuildNeighborTable();

protected:
	FEContactPotentialSurface	m_surf1;
	FEContactPotentialSurface	m_surf2;

protected:
	double	m_kc;
	double	m_p;
	double	m_Rin;
	double	m_Rout;
	double	m_Rmin;
	double	m_wtol;

	double	m_c1, m_c2;

	std::vector< std::set<FESurfaceElement*> >	m_activeElements;
	std::vector< std::set<FESurfaceElement*> >	m_elemNeighbors;

	DECLARE_FECORE_CLASS();
};

