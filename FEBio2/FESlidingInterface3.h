#pragma once
#include "FEContactInterface.h"
#include "FECore/FESurface.h"
#include "FECore/vec2d.h"

//-----------------------------------------------------------------------------
class FESlidingSurface3 : public FEContactSurface
{
public:
	//! constructor
	FESlidingSurface3(FEM* pfem);
	
	//! destructor
	~FESlidingSurface3() {}
	
	//! initialization
	void Init();
	
	//! shallow copy
	void ShallowCopy(FESlidingSurface3& s);
	
	//! calculate the nodal normals
	void UpdateNodeNormals();
	
	void Serialize(DumpFile& ar);
	
	void SetPoroMode(bool bporo) { m_bporo = bporo; }
	
protected:
	FEM*	m_pfem;
	
public:
	bool						m_bporo;	//!< set poro-mode
	bool						m_bsolu;	//!< set solute-mode
	
	vector<double>				m_gap;	//!< gap function at integration points
	vector<vec3d>				m_nu;	//!< normal at integration points
	vector<vec2d>				m_rs;	//!< natural coordinates of projection of integration point
	vector<double>				m_Lmd;	//!< lagrange multipliers for displacements
	vector<double>				m_Lmp;	//!< lagrange multipliers for fluid pressures
	vector<double>				m_Lmc;	//!< lagrange multipliers for solute concentrations
	vector<FESurfaceElement*>	m_pme;	//!< master element of projected integration point
	vector<int>					m_nei;	//!< surface element indices into arrays
	
	vector<double>	m_epsn;	//!< penalty factors
	vector<double>	m_epsp;	//!< pressure penalty factors
	vector<double>	m_epsc;	//!< concentration penalty factors
	
	vector<vec3d>		m_nn;	//!< node normals
	
	// biphasic-solute data
	vector<double>				m_pg;	//!< pressure "gap"
	vector<double>				m_cg;	//!< concentration "gap"
};

//-----------------------------------------------------------------------------
class FESlidingInterface3 :	public FEContactInterface
{
public:
	//! constructor
	FESlidingInterface3(FEM* pfem);
	
	//! destructor
	~FESlidingInterface3();
	
	//! initialization
	void Init();
	
	//! update
	void Update();
	
	//! Create a shallow copy
	void ShallowCopy(FEContactInterface& ci);
	
	//! calculate contact forces
	void ContactForces(vector<double>& F);
	
	//! calculate contact stiffness
	void ContactStiffness();
	
	//! calculate Lagrangian augmentations
	bool Augment(int naug);
	
	//! serialize data to archive
	void Serialize(DumpFile& ar);
	
	//! mark ambient condition 
	void MarkAmbient();
	
	//! set ambient condition 
	void SetAmbient();
	
protected:
	void ProjectSurface(FESlidingSurface3& ss, FESlidingSurface3& ms, bool bupseg);
	
	//! calculate penalty factor
	void CalcAutoPenalty(FESlidingSurface3& s);
	
	void CalcAutoPressurePenalty(FESlidingSurface3& s);
	double AutoPressurePenalty(FESurfaceElement& el, FESlidingSurface3& s);
	
	void CalcAutoConcentrationPenalty(FESlidingSurface3& s);
	double AutoConcentrationPenalty(FESurfaceElement& el, FESlidingSurface3& s);

	void BiphasicSoluteStatus(FEMesh& m, FESurfaceElement& el, bool& bstat, bool& sstat);
	
public:
	FESlidingSurface3	m_ms;	//!< master surface
	FESlidingSurface3	m_ss;	//!< slave surface
	
	int				m_knmult;	//!< higher order stiffness multiplier
	int				m_npass;	//!< nr of passes
	double			m_atol;		//!< augmentation tolerance
	double			m_gtol;		//!< gap tolerance
	double			m_ptol;		//!< pressure gap tolerance
	double			m_ctol;		//!< concentration gap tolerance
	double			m_stol;		//!< search tolerance
	bool			m_bsymm;	//!< use symmetric stiffness components only
	double			m_srad;		//!< contact search radius
	int				m_naugmax;	//!< maximum nr of augmentations
	int				m_naugmin;	//!< minimum nr of augmentations
	int				m_nsegup;	//!< segment update parameter
	
	double			m_epsn;		//!< normal penalty factor
	bool			m_bautopen;	//!< use autopenalty factor
	
	bool	m_bdebug;		// debug flag
	char	m_szdebug[256];	// debug file name
	FILE*	m_fp;			// debug file
	
	// biphasic-solute contact parameters
	double	m_epsp;		//!< fluid volumetric flow rate penalty
	double	m_epsc;		//!< solute molar flow rate penalty
	double	m_Rgas;		//!< universal gas constant
	double	m_Tabs;		//!< absolute temperature
	double	m_ambp;		//!< ambient pressure
	double	m_ambc;		//!< ambient concentration
	int		m_aplc;		//!< ambient pressure load curve
	int		m_aclc;		//!< ambient concentration load curve
	FELoadCurve* m_pplc;	//!< pointer to ambient pressure load curve
	FELoadCurve* m_pclc;	//!< pointer to ambient concentration load curve
};
