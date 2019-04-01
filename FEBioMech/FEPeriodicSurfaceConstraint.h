/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, Columbia University, and others.

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

//-----------------------------------------------------------------------------
class FEBIOMECH_API FEPeriodicSurfaceConstraintSurface : public FEContactSurface
{
public:
	//! constructor
	FEPeriodicSurfaceConstraintSurface(FEModel* pfem) : FEContactSurface(pfem) { m_nref = -1; }

	//! initializes data
	bool Init();

	//! calculates the center of mass of the surface
	vec3d CenterOfMass();

	void Serialize(DumpStream& ar);

public:
	vector<vec3d>				m_gap;	//!< gap function at nodes
	vector<FESurfaceElement*>	m_pme;	//!< master element a slave node penetrates
	vector<vec2d>				m_rs;	//!< natural coordinates of slave projection on master element
	vector<vec3d>				m_Lm;	//!< Lagrange multipliers

	int		m_nref;	//!< reference node
};

//-----------------------------------------------------------------------------

class FEBIOMECH_API FEPeriodicSurfaceConstraint : public FEContactInterface
{
public:
	//! constructor
	FEPeriodicSurfaceConstraint(FEModel* pfem);

	//! destructor
	virtual ~FEPeriodicSurfaceConstraint(void) {}

	//! initialization
	bool Init() override;

	//! interface activation
	void Activate() override;

	//! serialize data to archive
	void Serialize(DumpStream& ar) override;

	//! return the master and slave surface
	FESurface* GetMasterSurface() override { return &m_ms; }
	FESurface* GetSlaveSurface() override { return &m_ss; }

	//! return integration rule class
	bool UseNodalIntegration() override { return true; }

	//! build the matrix profile for use in the stiffness matrix
	void BuildMatrixProfile(FEGlobalMatrix& K) override;

public:
	//! calculate contact forces
	void Residual(FEGlobalVector& R, const FETimeInfo& tp) override;

	//! calculate contact stiffness
	void StiffnessMatrix(FESolver* psolver, const FETimeInfo& tp) override;

	//! calculate Lagrangian augmentations
	bool Augment(int naug, const FETimeInfo& tp) override;

	//! update
	void Update() override;

protected:
	void ProjectSurface(FEPeriodicSurfaceConstraintSurface& ss, FEPeriodicSurfaceConstraintSurface& ms, bool bmove);

public:
	FEPeriodicSurfaceConstraintSurface	m_ss;	//!< slave surface
	FEPeriodicSurfaceConstraintSurface	m_ms;	//!< master surface

	double	m_atol;			//!< augmentation tolerance
	double	m_eps;			//!< penalty scale factor
	double	m_stol;			//!< search tolerance
	double	m_srad;			//!< search radius (%)
	bool	m_btwo_pass;	//!< nr of passes

	int	m_dofX;
	int m_dofY;
	int	m_dofZ;

	DECLARE_FECORE_CLASS();
};
