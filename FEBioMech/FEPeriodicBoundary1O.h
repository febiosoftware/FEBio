#pragma once

#include "FEContactInterface.h"
#include "FEPeriodicBoundary.h"
#include "FEContactSurface.h"
#include <FECore/tens3d.h>

//-----------------------------------------------------------------------------

class FEPeriodicBoundary1O : public FEContactInterface
{
public:
	//! constructor
	FEPeriodicBoundary1O(FEModel* pfem);

	//! destructor
	virtual ~FEPeriodicBoundary1O(void) {}

	//! initialization
	bool Init();

	//! interface activation
	void Activate();

	//! update
	void Update(int niter);

	//! calculate contact forces
	void ContactForces(FEGlobalVector& R);

	//! calculate contact stiffness
	void ContactStiffness(FESolver* psolver);

	//! calculate Lagrangian augmentations
	bool Augment(int naug);

	//! serialize data to archive
	void Serialize(DumpStream& ar);

	//! return the master and slave surface
	FESurface* GetMasterSurface() { return &m_ms; }
	FESurface* GetSlaveSurface () { return &m_ss; }

	//! return integration rule class
	bool UseNodalIntegration() { return true; }

	//! build the matrix profile for use in the stiffness matrix
	void BuildMatrixProfile(FEGlobalMatrix& K);

	//! create a copy of this interface
	void CopyFrom(FESurfacePairInteraction* pci);

protected:
	void ProjectSurface(FEPeriodicSurface& ss, FEPeriodicSurface& ms, bool bmove);

public:
	FEPeriodicSurface		m_ss;	//!< slave surface
	FEPeriodicSurface		m_ms;	//!< master surface

	double	m_atol;			//!< augmentation tolerance
	double	m_eps;			//!< penalty scale factor
	double	m_stol;			//!< search tolerance
	double  m_srad;			//!< search radius (%)
	bool	m_btwo_pass;	//!< two-pass flag
	int		m_naugmin;		//!< minimum number of augmentations
	vec3d	m_off;			//!< relative displacement offset

	mat3d		m_Fmacro;		//!< Macroscopic deformation gradient
	
	DECLARE_PARAMETER_LIST();
};
