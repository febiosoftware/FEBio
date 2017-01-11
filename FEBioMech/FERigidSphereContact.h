#pragma once
#include <FECore/vector.h>
#include <FECore/FESurface.h>
#include <FECore/vec3d.h>
#include <FECore/vec2d.h>
#include <FECore/FENNQuery.h>
#include "FEContactInterface.h"
#include "FERigidSurface.h"

//-----------------------------------------------------------------------------
class FERigidSphereSurface : public FESurface
{
public:
	struct DATA
	{
		double	gap;	//!< gap function
		vec3d	nu;		//!< master normal at slave point
		double	Lm;		//!< Lagrange multiplier
	};

public:
	//! constructor
	FERigidSphereSurface(FEModel* pfem);

	//! Initializes data structures
	bool Init();

	//! Update the surface data
	void Update() {}

	//! Calculate the total traction at a node
	vec3d traction(int inode);

	// serialize surface data
	void Serialize(DumpStream& ar);

	//! Unpack surface element data
	void UnpackLM(FEElement& el, vector<int>& lm);

public:
	int	m_dofX;
	int	m_dofY;
	int	m_dofZ;

	vector<DATA>	m_data;		//!< integration point data
	FENNQuery		m_NQ;		//!< this structure is used in finding the master element that corresponds to a slave node
};

//-----------------------------------------------------------------------------
//! This class implements a sliding contact interface with a rigid wall

//! This class is a specialization of the general sliding interface where
//! the master surface is a rigid wall

class FERigidSphereContact : public FEContactInterface
{
public:
	//! constructor
	FERigidSphereContact(FEModel* pfem);

	//! intializes rigid wall interface
	bool Init();

	//! interface activation
	void Activate();

	//! project slave nodes onto master plane
	void ProjectSurface(FERigidSphereSurface& s);

	//! update rigid wall data
	void Update(int niter);

	//! calculate contact forces
	virtual void ContactForces(FEGlobalVector& R);

	//! calculate contact stiffness
	virtual void ContactStiffness(FESolver* psolver);

	//! calculate Lagrangian augmentations
	virtual bool Augment(int naug);

	//! serialize data to archive
	void Serialize(DumpStream& ar);

	//! return the master and slave surface
	FESurface* GetMasterSurface() { return 0; }
	FESurface* GetSlaveSurface() { return &m_ss; }

	//! return integration rule class
	bool UseNodalIntegration() { return true; }

	//! build the matrix profile for use in the stiffness matrix
	void BuildMatrixProfile(FEGlobalMatrix& K);

public:
	FERigidSphereSurface	m_ss;		//!< slave surface
	FERigidSphere			m_sphere;	//!< master surface

	int nse;	//!< number of slave elements

	double		m_atol;		//!< augmentation tolerance
	double		m_eps;		//!< penalty scale factor

	DECLARE_PARAMETER_LIST();
};
