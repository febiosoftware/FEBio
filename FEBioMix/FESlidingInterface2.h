#pragma once
#include "FEBioMech/FEContactInterface.h"
#include "FEBiphasicContactSurface.h"

//-----------------------------------------------------------------------------
class FESlidingSurface2 : public FEBiphasicContactSurface
{
public:
	//! Integration point data
	class Data : public FEBiphasicContactPoint
	{
	public:
		Data();

	public:
		double	m_Lmd;	//!< lagrange multipliers for displacement
		double	m_epsn;	//!< penalty factor
		double	m_epsp;	//!< pressure penatly factor
        double  m_p1;   //!< fluid pressure
		vec3d	m_nu;	//!< normal at integration points
		vec2d	m_rs;	//!< natrual coordinates of projection
		FESurfaceElement*	m_pme;	//!< master element
	};

public:
	//! constructor
	FESlidingSurface2(FEModel* pfem);

	//! initialization
	bool Init();

	// data serialization
	void Serialize(DumpStream& ar);

	//! evaluate net contact force
	vec3d GetContactForce();
    vec3d GetContactForceFromElementStress();

	//! evaluate net contact area
	double GetContactArea();
    
	//! evaluate net fluid force
	vec3d GetFluidForce();
    vec3d GetFluidForceFromElementPressure();
    
    //! evaluate the fluid load support
    double GetFluidLoadSupport();

	//! calculate the nodal normals
	void UpdateNodeNormals();

	void SetPoroMode(bool bporo) { m_bporo = bporo; }

	//! create material point data
	FEMaterialPoint* CreateMaterialPoint() override;

public:
    void GetContactTraction(int nface, vec3d& pt);
	void GetNodalContactPressure(int nface, double* pg);
	void GetNodalContactTraction(int nface, vec3d* pt);
    void EvaluateNodalContactPressures();

private:
	void GetContactPressure(int nface, double& pg);

protected:
	FEModel*	m_pfem;

public:
	bool	m_bporo;	//!< set poro-mode

	vector<bool>		m_poro;	//!< surface element poro status
	vector<vec3d>		m_nn;	//!< node normals
    vector<double>      m_pn;   //!< nodal contact pressures

	vec3d	m_Ft;	//!< total contact force (from equivalent nodal forces)
};

//-----------------------------------------------------------------------------
class FESlidingInterface2 :	public FEContactInterface
{
public:
	//! constructor
	FESlidingInterface2(FEModel* pfem);

	//! destructor
	~FESlidingInterface2();

	//! initialization
	bool Init() override;

	//! interface activation
	void Activate() override;

	//! calculate contact pressures for file output
	void UpdateContactPressures();

	//! serialize data to archive
	void Serialize(DumpStream& ar) override;

	//! mark free-draining condition 
	void MarkFreeDraining();
	
	//! set free-draining condition 
	void SetFreeDraining();

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
	void ProjectSurface(FESlidingSurface2& ss, FESlidingSurface2& ms, bool bupseg, bool bmove = false);

	//! calculate penalty factor
	void CalcAutoPenalty(FESlidingSurface2& s);

	void CalcAutoPressurePenalty(FESlidingSurface2& s);
	double AutoPressurePenalty(FESurfaceElement& el, FESlidingSurface2& s);

public:
	FESlidingSurface2	m_ms;	//!< master surface
	FESlidingSurface2	m_ss;	//!< slave surface

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
	int				m_nsegup;		//!< segment update parameter
	bool			m_breloc;		//!< node relocation on startup
    bool            m_bsmaug;       //!< smooth augmentation
    bool            m_bdupr;        //!< dual projection flag for free-draining

	double			m_epsn;		//!< normal penalty factor
	bool			m_bautopen;	//!< use autopenalty factor

	// biphasic contact parameters
	double	m_epsp;		//!< flow rate penalty

protected:
	int	m_dofP;

	DECLARE_FECORE_CLASS();
};
