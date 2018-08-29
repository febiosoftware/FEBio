#pragma once
#include "FEBioMech/FEContactInterface.h"
#include "FEBiphasicContactSurface.h"

//-----------------------------------------------------------------------------
class FETiedBiphasicSurface : public FEBiphasicContactSurface
{
public:
    //! Integration point data
    class Data
    {
    public:
        Data();
        
    public:
        vec3d	m_Gap;	//!< initial gap in reference configuration
        vec3d	m_dg;	//!< gap function at integration points
        vec3d	m_nu;	//!< normal at integration points
        vec2d	m_rs;	//!< natural coordinates of projection of integration point
        vec3d	m_Lmd;	//!< lagrange multipliers for displacements
        vec3d   m_tr;   //!< contact traction
        double	m_Lmp;	//!< lagrange multipliers for fluid pressures
        double	m_epsn;	//!< penalty factors
        double	m_epsp;	//!< pressure penalty factors
        double	m_pg;	//!< pressure "gap"
        FESurfaceElement*	m_pme;	//!< master element of projected integration point
    };
    
	//! constructor
	FETiedBiphasicSurface(FEModel* pfem);
	
	//! initialization
	bool Init();
	
	//! calculate the nodal normals
	void UpdateNodeNormals();
	
	void Serialize(DumpStream& ar);
	
	void SetPoroMode(bool bporo) { m_bporo = bporo; }
	
public:
    void GetVectorGap      (int nface, vec3d& pg);
    void GetContactTraction(int nface, vec3d& pt);
    
protected:
	FEModel*	m_pfem;
	
public:
	bool				m_bporo;	//!< set poro-mode
	
    vector< vector<Data> >  m_Data;	//!< integration point data
    vector<bool>			m_poro;	//!< surface element poro status
    vector<vec3d>			m_nn;	//!< node normals
};

//-----------------------------------------------------------------------------
class FETiedBiphasicInterface :	public FEContactInterface
{
public:
	//! constructor
	FETiedBiphasicInterface(FEModel* pfem);
	
	//! destructor
	~FETiedBiphasicInterface();
	
	//! initialization
	bool Init() override;
	
	//! interface activation
	void Activate() override;

	//! serialize data to archive
	void Serialize(DumpStream& ar) override;

	//! return the master and slave surface
	FESurface* GetMasterSurface() override { return &m_ms; }
	FESurface* GetSlaveSurface () override { return &m_ss; }

	//! return integration rule class
	bool UseNodalIntegration() override { return false; }

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
	void Update(int niter, const FETimeInfo& tp) override;

protected:
	void InitialProjection(FETiedBiphasicSurface& ss, FETiedBiphasicSurface& ms);
	void ProjectSurface(FETiedBiphasicSurface& ss, FETiedBiphasicSurface& ms);
	
	//! calculate penalty factor
	void CalcAutoPenalty(FETiedBiphasicSurface& s);
	
	void CalcAutoPressurePenalty(FETiedBiphasicSurface& s);
	double AutoPressurePenalty(FESurfaceElement& el, FETiedBiphasicSurface& s);
	
public:
	FETiedBiphasicSurface	m_ms;	//!< master surface
	FETiedBiphasicSurface	m_ss;	//!< slave surface
	
	int				m_knmult;		//!< higher order stiffness multiplier
	bool			m_btwo_pass;	//!< two-pass flag
	double			m_atol;			//!< augmentation tolerance
	double			m_gtol;			//!< gap tolerance
	double			m_ptol;			//!< pressure gap tolerance
	double			m_stol;			//!< search tolerance
	bool			m_bsymm;		//!< use symmetric stiffness components only
	double			m_srad;			//!< contact search radius
	int				m_naugmax;		//!< maximum nr of augmentations
	int				m_naugmin;		//!< minimum nr of augmentations
	
	double			m_epsn;			//!< normal penalty factor
	bool			m_bautopen;		//!< use autopenalty factor
	
	// biphasic contact parameters
	double			m_epsp;		//!< flow rate penalty

protected:
	int	m_dofP;
	
	DECLARE_PARAMETER_LIST();
};
