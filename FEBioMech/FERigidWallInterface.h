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
#include <FECore/FESurfaceConstraint.h>
#include <FECore/FESurface.h>
#include <FECore/vec3d.h>
#include <FECore/vec2d.h>
#include <FECore/FENNQuery.h>

//-----------------------------------------------------------------------------
class FERigidWallSurface : public FESurface
{
public:
	//! constructor
	FERigidWallSurface(FEModel* pfem);

	//! Initializes data structures
	bool Init();

	//! Update the surface data
	void Update() {}

	//! Calculate the total traction at a node
	vec3d traction(int inode);

	void UpdateNormals();

	void Serialize(DumpStream& ar);

	//! Unpack surface element data
	void UnpackLM(FEElement& el, vector<int>& lm);

public:
	vector<double>				m_gap;	//!< gap function at nodes
	vector<vec3d>				m_nu;	//!< secondary surface normal at primary surface node
	vector<FESurfaceElement*>	m_pme;	//!< secondary surface element a primary surface node penetrates
	vector<vec2d>				m_rs;	//!< natural coordinates of projection on secondary surface element
	vector<vec2d>				m_rsp;	//!< natural coordinates at previous time step
	vector<double>				m_Lm;	//!< Lagrange multipliers for contact pressure
	vector<mat2d>				m_M;	//!< surface metric tensor
	vector<vec2d>				m_Lt;	//!< Lagrange multipliers for friction
	vector<double>				m_off;	//!< gap offset (= shell thickness)
	vector<double>				m_eps;	//!< normal penalty factors

	int	m_dofX;
	int	m_dofY;
	int	m_dofZ;

	FENNQuery		m_NQ;		//!< used in finding the secondary surface element that corresponds to a primary surface node
};

//-----------------------------------------------------------------------------
//! This class implements a sliding contact interface with a rigid wall

//! This class is a specialization of the general sliding interface where
//! the secondary surface is a rigid wall

class FERigidWallInterface : public FESurfaceConstraint
{
public:
	//! constructor
	FERigidWallInterface(FEModel* pfem);

	~FERigidWallInterface();

	//! intializes rigid wall interface
	bool Init() override;

	//! interface activation
	void Activate() override;

	//! serialize data to archive
	void Serialize(DumpStream& ar) override;

	//! return the primary and secondary surface
	FESurface* GetSurface() override { return m_ss; }

	//! return integration rule class
	//! TODO: Need to fix this! Maybe move to FESurfaceConstraint?
//	bool UseNodalIntegration() override { return true; }

	//! build the matrix profile for use in the stiffness matrix
	void BuildMatrixProfile(FEGlobalMatrix& K) override;

public:
	//! calculate contact forces
	void LoadVector(FEGlobalVector& R, const FETimeInfo& tp) override;

	//! calculate contact stiffness
	void StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp) override;

	//! calculate Lagrangian augmentations
	bool Augment(int naug, const FETimeInfo& tp) override;

	//! update
	void Update() override;

private:
	//! project surface nodes onto plane
	void ProjectSurface(FERigidWallSurface& s);

private:
	//! return plane normal
	vec3d PlaneNormal(const vec3d& r);

	//! project node onto plane
	vec3d ProjectToPlane(const vec3d& r);

public:
	FERigidWallSurface*	m_ss;		//!< primary surface

	double	m_a[4];	//!< plane equation

	int			m_laugon;	//!< augmentation flag
	double		m_atol;		//!< augmentation tolerance
	double		m_eps;		//!< penalty scale factor
	double		m_d;		//!< normal offset

	DECLARE_FECORE_CLASS();
};
