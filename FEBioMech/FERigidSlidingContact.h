#pragma once
#include <FECore/vector.h>
#include <FECore/FESurface.h>
#include <FECore/vec3d.h>
#include <FECore/vec2d.h>
#include <FECore/FENNQuery.h>
#include "FEContactInterface.h"
#include "FEContactSurface.h"
#include "FERigidSurface.h"

//-----------------------------------------------------------------------------
class FERigidSlidingSurface : public FEContactSurface
{
public:
	struct DATA
	{
		double	gap;	//!< gap function
		vec3d	nu;		//!< master normal at slave point
		double	Lm;		//!< Lagrange multiplier
		double	eps;	//!< penalty parameter
	};

public:
	//! constructor
	FERigidSlidingSurface(FEModel* pfem);

	//! Initializes data structures
	bool Init();

	//! Update the surface data
	void Update() {}

	// serialize surface data
	void Serialize(DumpStream& ar);

	//! Unpack surface element data
	void UnpackLM(FEElement& el, vector<int>& lm);

	//! Calculate the total traction at a node
	vec3d traction(int inode);

	//! total contact force
	vec3d GetContactForce();

public:
	int	m_dofX;
	int	m_dofY;
	int	m_dofZ;
	vec3d	m_Fc;	//!< total contact force

	vector<DATA>	m_data;		//!< integration point data
	FENNQuery		m_NQ;		//!< this structure is used in finding the master element that corresponds to a slave node
};

//-----------------------------------------------------------------------------
//! This class implements a sliding contact interface with a rigid boundary.
//! The rigid boundary is represented by a rigid surface
//
class FERigidSlidingContact : public FEContactInterface
{
public:
	//! constructor
	FERigidSlidingContact(FEModel* pfem);

	//! destructor
	~FERigidSlidingContact();

	//! intializes rigid sliding interface
	bool Init() override;

	//! interface activation
	void Activate() override;

	//! project slave nodes onto master plane
	void ProjectSurface(FERigidSlidingSurface& s);

	//! serialize data to archive
	void Serialize(DumpStream& ar) override;

	//! return the master and slave surface
	FESurface* GetMasterSurface() override { return 0; }
	FESurface* GetSlaveSurface() override { return &m_ss; }

	//! return integration rule class
	bool UseNodalIntegration() override { return false; }

	//! build the matrix profile for use in the stiffness matrix
	void BuildMatrixProfile(FEGlobalMatrix& K) override;

	//! Calculate auto-penalty
	void CalcAutoPenalty(FERigidSlidingSurface& s);

public:
	//! calculate contact forces
	void Residual(FEGlobalVector& R, const FETimeInfo& tp) override;

	//! calculate contact stiffness
	void StiffnessMatrix(FESolver* psolver, const FETimeInfo& tp) override;

	//! calculate Lagrangian augmentations
	bool Augment(int naug, const FETimeInfo& tp) override;

	//! update
	void Update(int niter, const FETimeInfo& tp) override;

private:
	FERigidSlidingSurface	m_ss;		//!< slave surface
	FERigidSurface*			m_rigid;	//!< master surface


public: // parameters
	double		m_atol;				//!< augmentation tolerance
	double		m_eps;				//!< penalty scale factor
	double		m_gtol;				//!< gap tolerance
	int			m_naugmin;			//!< min nr of augmentations
	int			m_naugmax;			//!< max nr of augmentations
	bool		m_bautopen;			//!< auto-penalty flag
	std::string	m_rigidName;		//!< name of rigid surface object

	DECLARE_PARAMETER_LIST();
};
