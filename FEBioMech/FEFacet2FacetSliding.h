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
//! Contact surface for facet-to-facet sliding interfaces
class FEFacetSlidingSurface : public FEContactSurface
{
public:
	//! Integration point data
	class Data : public FEContactMaterialPoint
	{
	public:
		Data();

		void Serialize(DumpStream& ar);

	public:
		double	m_Lm;	//!< Lagrange multipliers
		double	m_eps;	//!< penalty value at integration point
		vec3d	m_nu;	//!< secondary surface normal at integration points
		vec2d	m_rs;	//!< natural coordinates of projection of integration point
	};

public:
	//! constructor
	FEFacetSlidingSurface(FEModel* pfem);

	//! initialization
	bool Init() override;

	//! evaluate net contact force
	vec3d GetContactForce() override;

	//! evaluate net contact area
	double GetContactArea() override;
    
	//! create material point data
	FEMaterialPoint* CreateMaterialPoint() override;

	//! serialization
	void Serialize(DumpStream& ar) override;

public:
    void GetContactTraction(int nface, vec3d& pt) override;
	void GetNodalContactPressure(int nface, double* pn) override;
	void GetNodalContactTraction(int nface, vec3d* tn) override;

public:
	vector<vec3d>	m_Fn;	//!< equivalent nodal forces
};

//-----------------------------------------------------------------------------
//! Sliding interface with facet to facet integration

//! This class is similar to the sliding interface except that it uses
//! a Gaussian quadrature rule in stead of a nodal integration rule
//! as its base class does.
//
class FEFacet2FacetSliding : public FEContactInterface
{
public:
	//! constructor
	FEFacet2FacetSliding(FEModel* pfem);

	//! initialization routine
	bool Init() override;

	//! interface activation
	void Activate() override;

	//! calculate contact pressures for file output
	void UpdateContactPressures();

	//! serialize data to archive
	void Serialize(DumpStream& ar) override;

	//! get primary and secondary surfaces
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

	//! update 
	void Update() override;

protected:
	//! project primary surface onto secondary
	void ProjectSurface(FEFacetSlidingSurface& ss, FEFacetSlidingSurface& ms, bool bsegup, bool bmove = false);

	//! calculate auto-penalty
    void UpdateAutoPenalty();
    
	void CalcAutoPenalty(FEFacetSlidingSurface& s);

public:
	double	m_epsn;			//!< normal penalty factor
	double	m_knmult;		//!< normal stiffness multiplier
	double	m_stol;			//!< search tolerance
	bool	m_btwo_pass;	//!< two-pass flag
	bool	m_bautopen;		//!< auto-penalty flag
    bool    m_bupdtpen;     //!< update penalty at each time step
	double	m_srad;			//!< search radius (% of model size)
	int		m_nsegup;		//!< segment update parameter
	bool	m_breloc;       //!< node relocation on initialization
    bool    m_bsmaug;       //!< smooth augmentation

	double	m_atol;			//!< aug lag tolerance
	double	m_gtol;			//!< gap tolerance
	int		m_naugmin;		//!< min nr of augmentations
	int		m_naugmax;		//!< max nr of augmentations

	double			m_mu;		//!< friction coefficient (not implemented yet)
	double			m_epsf;		//!< penalty scale factor for friction (not implementer yet)

	double	m_dxtol;		//!< penalty insertion distance

	FEFacetSlidingSurface	m_ss;	//!< primary surface
	FEFacetSlidingSurface	m_ms;	//!< secondary surface

private:
	bool	m_bfirst;
	double	m_normg0;

public:
	DECLARE_FECORE_CLASS();
};
