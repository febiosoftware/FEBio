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
#include "FESolute.h"
#include <FECore/FECoreClass.h>
#include <map>

//-----------------------------------------------------------------------------
class FEBIOMIX_API FESlidingSurfaceMP : public FEBiphasicContactSurface
{
public:
    //! integration point data
    class Data : public FEBiphasicContactPoint
    {
    public:
        Data();
        
        void Init() override;
        
        void Serialize(DumpStream& ar) override;
        
    public:
        vec3d   m_dg;       //!< vector gap
        double  m_Lmd;      //!< Lagrange multipliers for displacements
        vec3d   m_Lmt;      //!< Lagrange multipliers for vector traction
        double  m_epsn;     //!< displacement penalty factors
        double  m_epsp;     //!< pressure penalty factors
        double  m_p1;       //!< fluid pressure
        vec3d   m_nu;       //!< normal at integration points
        vec3d   m_s1;       //!< tangent along slip direction
        vec3d   m_tr;       //!< contact traction
        vec2d   m_rs;       //!< natural coordinates of projection of integration point
        vec2d   m_rsp;      //!< m_rs at the previous time step
        bool    m_bstick;   //!< stick flag
        vector<double>  m_Lmc;  //!< Lagrange multipliers for solute concentrations
        vector<double>  m_epsc; //!< concentration penalty factors
        vector<double>  m_cg;   //!< concentration "gap"
        vector<double>  m_c1;   //!< solute concentration
    };

public:
	//! constructor
	FESlidingSurfaceMP(FEModel* pfem);
	
	//! destructor
	~FESlidingSurfaceMP();
	
	//! initialization
	bool Init() override;
	
    //! initialize sliding surface and store previous values
    void InitSlidingSurface();
    
	//! evaluate net contact force
	vec3d GetContactForce() override;
	
	//! evaluate net contact area
	double GetContactArea() override;
    
	//! evaluate net fluid force
	vec3d GetFluidForce() override;
	
	//! calculate the nodal normals
	void UpdateNodeNormals();
	
	void Serialize(DumpStream& ar) override;
	
	void SetPoroMode(bool bporo) { m_bporo = bporo; }

	void UnpackLM(FEElement& el, vector<int>& lm) override;
	
	//! create material point data
	FEMaterialPoint* CreateMaterialPoint() override;

public:
    void GetVectorGap           (int nface, vec3d& pg) override;
    void GetContactTraction     (int nface, vec3d& pt) override;
    void GetSlipTangent         (int nface, vec3d& pt);
    void GetMuEffective         (int nface, double& pg) override;
    void GetLocalFLS            (int nface, double& pg) override;
    void GetNodalVectorGap      (int nface, vec3d* pg) override;
    void GetNodalContactPressure(int nface, double* pg) override;
    void GetNodalContactTraction(int nface, vec3d* tn) override;
    void GetStickStatus         (int nface, double& pg) override;
    void EvaluateNodalContactPressures();
    void EvaluateNodalContactTractions();

private:
	void GetContactPressure(int nface, double& pg);
	
public:
	bool						m_bporo;	//!< set poro-mode
	bool						m_bsolu;	//!< set solute-mode
	
	vector<vec3d>				m_nn;	//!< node normals
    vector<double>              m_pn;   //!< nodal contact pressures
    vector<vec3d>               m_tn;   //!< nodal contact tractions

	vector<int>					m_sid;	//!< list of solute id's for this surface
    
    vec3d	m_Ft;                       //!< total contact force (from equivalent nodal forces)

protected:
	int	m_dofC;
};

//-----------------------------------------------------------------------------
// helper class for reading ambient concentrations
class FEAmbientConcentration : public FECoreClass
{
public:
	FEAmbientConcentration(FEModel* fem);

public:
	int			m_sol;
	double		m_ambc;

	DECLARE_FECORE_CLASS();
	FECORE_BASE_CLASS(FEAmbientConcentration);
};

//-----------------------------------------------------------------------------
class FEBIOMIX_API FESlidingInterfaceMP : public FEContactInterface
{

public:
	//! constructor
	FESlidingInterfaceMP(FEModel* pfem);
	
	//! destructor
	~FESlidingInterfaceMP();
	
	//! initialization
	bool Init() override;
	
	//! interface activation
	void Activate() override;

    //! calculate the slip direction on the primary surface
    vec3d SlipTangent(FESlidingSurfaceMP& ss, const int nel, const int nint, FESlidingSurfaceMP& ms, double& dh, vec3d& r);

    //! calculate contact traction
    vec3d ContactTraction(FESlidingSurfaceMP& ss, const int nel, const int n, FESlidingSurfaceMP& ms, double& pn);
    
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
 
	//! get solute data
    FESoluteData* FindSoluteData(int nid);

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
	void ProjectSurface(FESlidingSurfaceMP& ss, FESlidingSurfaceMP& ms, bool bupseg, bool bmove = false);
	
	//! calculate penalty factor
    void UpdateAutoPenalty();
    
	void CalcAutoPenalty(FESlidingSurfaceMP& s);
	
	void CalcAutoPressurePenalty(FESlidingSurfaceMP& s);
	double AutoPressurePenalty(FESurfaceElement& el, FESlidingSurfaceMP& s);
	
	void CalcAutoConcentrationPenalty(FESlidingSurfaceMP& s, const int isol);
	double AutoConcentrationPenalty(FESurfaceElement& el, FESlidingSurfaceMP& s, const int isol);
    
    double AutoPenalty(FESurfaceElement& el, FESurface &s);
	
public:
	FESlidingSurfaceMP	m_ss;	//!< primary surface
	FESlidingSurfaceMP	m_ms;	//!< secondary surface
	
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
    bool            m_bsmfls;       //!< smooth local fluid load support

	double			m_epsn;         //!< normal penalty factor
	bool			m_bautopen;     //!< use autopenalty factor
    bool            m_bupdtpen;     //!< update penalty at each time step
    
    double          m_mu;           //!< friction coefficient
    bool            m_bfreeze;      //!< freeze stick/slip status

	// multiphasic contact parameters
    double  m_phi;                  //!< solid-solid contact fraction
	double	m_epsp;					//!< fluid volumetric flow rate penalty
	double	m_epsc;					//!< solute molar flow rate penalty
	double	m_Rgas;					//!< universal gas constant
	double	m_Tabs;					//!< absolute temperature
	double	m_ambp;					//!< ambient pressure
	vector<double>	m_ambc;         //!< ambient concentration

	vector<FEAmbientConcentration*>	m_ambctmp;
	vector<int> m_sid;				//!< list of solute ids common to both contact surfaces
	vector<int> m_ssl;				//!< list of primary surface solutes common to both contact surfaces
	vector<int> m_msl;				//!< list of secondary surface solutes common to both contact surfaces
    vector<int> m_sz;               //!< charge number of solutes common to both contact surfaces

protected:
	int	m_dofP;
	int	m_dofC;
	
	DECLARE_FECORE_CLASS();
};
