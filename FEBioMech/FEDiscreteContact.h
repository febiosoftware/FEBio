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
#include "FEContactSurface.h"
#include "FEDeformableSpringDomain.h"

//-----------------------------------------------------------------------------
class FEDiscreteSet;

//-----------------------------------------------------------------------------
class FEBIOMECH_API FEDiscreteContactSurface : public FEContactSurface
{
public:
	//! constructor
	FEDiscreteContactSurface(FEModel* fem);

	//! Initialization
	bool Init();
};

//-----------------------------------------------------------------------------
class FEBIOMECH_API FEDiscreteContact : public FESurfaceConstraint
{
	struct NODE
	{
		int					nid;	//!< (local) node ID
		FESurfaceElement*	pe;		//!< secondary surface element
		double				gap;	//!< gap distance
		double				Lm;		//!< Lagrange multiplier
		vec3d				nu;		//!< normal at secondary surface projection
		vec3d				q;		//!< projection point
		double				proj[2];	//!< iso-parametric coordinates of projection point
	};

public:
	FEDiscreteContact(FEModel* pfem);

public:
	bool Init() override;
	void Activate() override;

	void LoadVector(FEGlobalVector& R, const FETimeInfo& tp) override;
	void StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp) override;
	bool Augment(int naug, const FETimeInfo& tp) override;
	void BuildMatrixProfile(FEGlobalMatrix& M) override;
	void Update(const FETimeInfo& tp);

	void SetDiscreteSet(FEDiscreteSet* pset);

	FESurface* GetSurface() override { return &m_surf; }

protected:
	void ProjectSurface(bool bupseg);
	void ContactNodalForce    (NODE& nodeData, FESurfaceElement& mel, vector<double>& fe);
	void ContactNodalStiffness(NODE& nodeData, FESurfaceElement& mel, matrix& ke);

protected:
	FEDiscreteContactSurface	m_surf;
	vector<NODE>	m_Node;
	double	m_normg0;
	bool	m_bfirst;

protected:
	bool	m_blaugon;	//!< augmentation flag
	double	m_altol;	//!< augmentation tolerance
	double	m_penalty;	//!< penalty parameter
	double	m_gaptol;	//!< gap tolerance
	int		m_naugmin;	//!< minimum number of augmentations
	int		m_naugmax;	//!< maximum number of augmentations
	int		m_nsegup;	//!< number of segment updates (or zero)

	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
class FEBIOMECH_API FEDiscreteContact2 : public FESurfaceConstraint
{
	struct NODE
	{
		int		node;				// node index (local ID into discrete domain)
		FESurfaceElement*	pe;		// secondary surface element
		double	proj[2];			// natural coordinates of projection
		vec3d	nu;					// normal on secondary surface
		vec3d	q;					// new position
	};

public:
	FEDiscreteContact2(FEModel* fem);

	bool Init() override;
	void Activate() override;

	void LoadVector(FEGlobalVector& R, const FETimeInfo& tp) override;
	void StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp) override;
	void BuildMatrixProfile(FEGlobalMatrix& M) override;
	void Update(const FETimeInfo& tp);
	bool Augment(int naug, const FETimeInfo& tp) override { return true; }

	void SetDiscreteDomain(FEDeformableSpringDomain2* dom) { m_dom = dom; }
	FESurface* GetSurface() override { return &m_surf; }

protected:
	void ProjectNodes();

protected:
	FEDiscreteContactSurface	m_surf;
	FEDeformableSpringDomain2*	m_dom;
	vector<NODE>	m_nodeData;

	DECLARE_FECORE_CLASS();
};
