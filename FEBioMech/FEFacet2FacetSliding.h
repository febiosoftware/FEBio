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

	public:
		double	m_Lm;	//!< Lagrange multipliers
		double	m_eps;	//!< penalty value at integration point
		vec3d	m_nu;	//!< master normal at integration points
		vec2d	m_rs;	//!< natural coordinates of projection of integration point
		FESurfaceElement*	m_pme;	//!< master element of projection
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
    
	//! serialize data for (cold) restart
	void Serialize(DumpStream& ar) override;

	//! create material point data
	FEMaterialPoint* CreateMaterialPoint() override;

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
	//! project slave surface onto master
	void ProjectSurface(FEFacetSlidingSurface& ss, FEFacetSlidingSurface& ms, bool bsegup, bool bmove = false);

	//! calculate auto-penalty
	void CalcAutoPenalty(FEFacetSlidingSurface& s);

public:
	double	m_epsn;			//!< normal penalty factor
	double	m_knmult;		//!< normal stiffness multiplier
	double	m_stol;			//!< search tolerance
	bool	m_btwo_pass;	//!< two-pass flag
	bool	m_bautopen;		//!< auto-penalty flag
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

	FEFacetSlidingSurface	m_ms;	//!< master surface
	FEFacetSlidingSurface	m_ss;	//!< slave surface

private:
	bool	m_bfirst;
	double	m_normg0;

public:
	DECLARE_FECORE_CLASS();
};
