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

//-----------------------------------------------------------------------------
//! Surface definition for the facet-to-facet tied interface
class FEFacetTiedSurface : public FEContactSurface
{
public:
	//! integration point data
	class Data : public FEContactMaterialPoint
	{
	public:
		Data();

		void Serialize(DumpStream& ar);

	public:
		vec3d	m_vgap;		//!< gap function
		vec3d	m_vgap0;	//!< initial gap function
		vec3d	m_Lm;	//!< Lagrange multiplier
		vec2d	m_rs;	//!< natural coordinates on secondary surface element
	};

public:
	//! constructor
	FEFacetTiedSurface(FEModel* pfem);

	//! Initialization
	bool Init() override;

	//! serialization for cold restarts
	void Serialize(DumpStream& ar) override;

	//! create material point data
	FEMaterialPoint* CreateMaterialPoint() override;
};

//-----------------------------------------------------------------------------
//! Tied contact interface with facet-to-facet integration
class FEFacet2FacetTied : public FEContactInterface
{
public:
	//! constructor
	FEFacet2FacetTied(FEModel* pfem);

	//! Initialization
	bool Init() override;

	//! interface activation
	void Activate() override;

	//! serialize data to archive
	void Serialize(DumpStream& ar) override;

	//! return the primary and secondary surface
	FESurface* GetPrimarySurface() override { return &m_ss; }
	FESurface* GetSecondarySurface() override { return &m_ms; }

	//! return integration rule class
	bool UseNodalIntegration() override { return false; }

	//! build the matrix profile for use in the stiffness matrix
	void BuildMatrixProfile(FEGlobalMatrix& K) override;

public:
	//! calculate contact forces
	void LoadVector(FEGlobalVector& R, const FETimeInfo& tp) override;

	//! calculate contact stiffness
	void StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp) override;

	//! calculate Lagrangian augmentations
	bool Augment(int naug, const FETimeInfo& tp) override;

	//! update contact data
	void Update() override;

protected:

	//! projects primary nodes onto secondary nodes
	void ProjectSurface(FEFacetTiedSurface& ss, FEFacetTiedSurface& ms);

private:
	FEFacetTiedSurface	m_ms;	//!< secondary surface
	FEFacetTiedSurface	m_ss;	//!< primary surface

public:
	double		m_atol;		//!< augmentation tolerance
	double		m_eps;		//!< penalty scale factor
	double		m_stol;		//!< search tolerance
	int			m_naugmax;	//!< maximum nr of augmentations
	int			m_naugmin;	//!< minimum nr of augmentations
	bool		m_gapoff;	//!< retain initial gap as offset

	DECLARE_FECORE_CLASS();
};
