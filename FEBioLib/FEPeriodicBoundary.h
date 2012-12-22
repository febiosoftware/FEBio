#pragma once

#include "FECore/FEContactInterface.h"
#include "FEContactSurface.h"

//-----------------------------------------------------------------------------
class FEPeriodicSurface : public FEContactSurface
{
public:
	//! constructor
	FEPeriodicSurface(FEMesh* pm = 0) : FEContactSurface(pm) {}

	//! initializes data
	bool Init();

	//! shallow copy
	void ShallowCopy(FEPeriodicSurface& s)
	{
		m_Lm = s.m_Lm;
		m_gap = s.m_gap;
	}

	//! update surface data
	void Update();

	//! calculates the center of mass of the surface
	vec3d CenterOfMass();

	void Serialize(DumpFile& ar);

public:
	vector<vec3d>				m_gap;	//!< gap function at nodes
	vector<FESurfaceElement*>	m_pme;	//!< master element a slave node penetrates
	vector<vec2d>				m_rs;	//!< natural coordinates of slave projection on master element
	vector<vec3d>				m_Lm;	//!< Lagrange multipliers
};

//-----------------------------------------------------------------------------

class FEPeriodicBoundary : public FEContactInterface
{
public:
	//! constructor
	FEPeriodicBoundary(FEModel* pfem);

	//! destructor
	virtual ~FEPeriodicBoundary(void) {}

	//! initialization
	bool Init();

	//! interface activation
	void Activate();

	//! update
	void Update(int niter);

	//! shallow copy
	void ShallowCopy(FEContactInterface& ci);

	//! calculate contact forces
	void ContactForces(FEGlobalVector& R);

	//! calculate contact stiffness
	void ContactStiffness(FENLSolver* psolver);

	//! calculate Lagrangian augmentations
	bool Augment(int naug);

	//! serialize data to archive
	void Serialize(DumpFile& ar);

	//! return the master and slave surface
	FESurface* GetMasterSurface() { return &m_ms; }
	FESurface* GetSlaveSurface () { return &m_ss; }

protected:
	void ProjectSurface(FEPeriodicSurface& ss, FEPeriodicSurface& ms, bool bmove);

public:
	FEPeriodicSurface	m_ss;	//!< slave surface
	FEPeriodicSurface	m_ms;	//!< master surface

	double	m_atol;			//!< augmentation tolerance
	double	m_eps;			//!< penalty scale factor
	double	m_stol;			//!< search tolerance
	double  m_srad;			//!< search radius (%)
	bool	m_btwo_pass;	//!< two-pass flag

	DECLARE_PARAMETER_LIST();
};
