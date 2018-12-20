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
	bool UseNodalIntegration() override { return true; }

	//! build the matrix profile for use in the stiffness matrix
	void BuildMatrixProfile(FEGlobalMatrix& K) override;

	//! create a copy of this interface
	void CopyFrom(FESurfacePairConstraint* pci) override;

public:
	//! calculate contact forces
	void Residual(FEGlobalVector& R, const FETimeInfo& tp) override;

	//! calculate contact stiffness
	void StiffnessMatrix(FESolver* psolver, const FETimeInfo& tp) override;

	//! calculate Lagrangian augmentations
	bool Augment(int naug, const FETimeInfo& tp) override;

	//! update
	void Update() override;

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
	
	DECLARE_FECORE_CLASS();
};
