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
#include "FEBioMech/FEContactInterface.h"
#include "FEBiphasicContactSurface.h"

//-----------------------------------------------------------------------------
class FEBIOMIX_API FESlidingSurface3 : public FEBiphasicContactSurface
{
public:
	//! integration point data
	class Data : public FEBiphasicContactPoint
	{
	public:
		Data();

		void Serialize(DumpStream& ar) override;

	public:
		double	m_Lmd;	//!< Lagrange multipliers for displacements
		double	m_Lmc;	//!< Lagrange multipliers for solute concentrations
		double	m_epsn;	//!< displacement penalty factors
		double	m_epsp;	//!< pressure penalty factors
		double	m_epsc;	//!< concentration penatly factors
		double	m_cg;	//!< concentration "gap"
		vec3d	m_nu;	//!< normal at integration points
		vec2d	m_rs;	//!< natural coordinates of projection of integration point
	};

public:
	//! constructor
	FESlidingSurface3(FEModel* pfem);
	
	//! destructor
	~FESlidingSurface3();
	
	//! initialization
	bool Init() override;
	
	//! evaluate net contact force
	vec3d GetContactForce() override;
	
	//! evaluate net contact area
	double GetContactArea() override;
    
	//! evaluate net fluid force
	vec3d GetFluidForce() override;
	
	//! calculate the nodal normals
	void UpdateNodeNormals();
	
	// data serialization
	void Serialize(DumpStream& ar) override;
	
	void SetPoroMode(bool bporo) { m_bporo = bporo; }

	void UnpackLM(FEElement& el, vector<int>& lm) override;

	//! create material point data
	FEMaterialPoint* CreateMaterialPoint() override;

public:
    void GetContactTraction(int nface, vec3d& pt) override;
	void GetNodalContactPressure(int nface, double* pg) override;
	void GetNodalContactTraction(int nface, vec3d* tn) override;
    void EvaluateNodalContactPressures();

private:
	void GetContactPressure(int nface, double& pg);
	
public:
	bool						m_bporo;	//!< set poro-mode
	bool						m_bsolu;	//!< set solute-mode
	
	vector<bool>				m_poro;	//!< surface element poro status
	vector<int>					m_solu;	//!< surface element solute id

	vector<vec3d>		m_nn;	//!< node normals
    vector<double>      m_pn;   //!< nodal contact pressures
    
    vec3d	m_Ft;	//!< total contact force (from equivalent nodal forces)

protected:
	int	m_dofC;
};

//-----------------------------------------------------------------------------
class FEBIOMIX_API FESlidingInterface3 :	public FEContactInterface
{
public:
	//! constructor
	FESlidingInterface3(FEModel* pfem);
	
	//! destructor
	~FESlidingInterface3();
	
	//! initialization
	bool Init() override;
	
	//! interface activation
	void Activate() override;

	//! calculate contact pressures for file output
	void UpdateContactPressures();

	//! serialize data to archive
	void Serialize(DumpStream& ar) override;
	
	//! mark ambient condition 
	void MarkAmbient();
	
	//! set ambient condition 
	void SetAmbient();

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

	//! update
	void Update() override;

protected:
	void ProjectSurface(FESlidingSurface3& ss, FESlidingSurface3& ms, bool bupseg, bool bmove = false);
	
	//! calculate penalty factor
    void UpdateAutoPenalty();
    
	void CalcAutoPenalty(FESlidingSurface3& s);
	
	void CalcAutoPressurePenalty(FESlidingSurface3& s);
	double AutoPressurePenalty(FESurfaceElement& el, FESlidingSurface3& s);
	
	void CalcAutoConcentrationPenalty(FESlidingSurface3& s);
	double AutoConcentrationPenalty(FESurfaceElement& el, FESlidingSurface3& s);

public:
	FESlidingSurface3	m_ss;	//!< primary surface
	FESlidingSurface3	m_ms;	//!< secondary surface
	
	int				m_knmult;		//!< higher order stiffness multiplier
	bool			m_btwo_pass;	//!< two-pass flag
	double			m_atol;			//!< augmentation tolerance
	double			m_gtol;			//!< gap tolerance
	double			m_ptol;			//!< pressure gap tolerance
	double			m_ctol;			//!< concentration gap tolerance
	double			m_stol;			//!< search tolerance
	bool			m_bsymm;		//!< use symmetric stiffness components only
	double			m_srad;			//!< contact search radius
	int				m_naugmax;		//!< maximum nr of augmentations
	int				m_naugmin;		//!< minimum nr of augmentations
	int				m_nsegup;		//!< segment update parameter
    bool			m_breloc;		//!< node relocation on startup
    bool            m_bsmaug;       //!< smooth augmentation
	
	double			m_epsn;		//!< normal penalty factor
	bool			m_bautopen;	//!< use autopenalty factor
    bool            m_bupdtpen;     //!< update penalty at each time step

	// biphasic-solute contact parameters
	double	m_epsp;		//!< fluid volumetric flow rate penalty
	double	m_epsc;		//!< solute molar flow rate penalty
	double	m_Rgas;		//!< universal gas constant
	double	m_Tabs;		//!< absolute temperature
	double	m_ambp;		//!< ambient pressure
	double	m_ambc;		//!< ambient concentration

protected:
	int	m_dofP;
	int	m_dofC;

	DECLARE_FECORE_CLASS();
};
