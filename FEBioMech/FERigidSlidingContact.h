/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2020 University of Utah, The Trustees of Columbia University in 
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

		void Serialize(DumpStream& ar);
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
	void LoadVector(FEGlobalVector& R, const FETimeInfo& tp) override;

	//! calculate contact stiffness
	void StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp) override;

	//! calculate Lagrangian augmentations
	bool Augment(int naug, const FETimeInfo& tp) override;

	//! update
	void Update() override;

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

	DECLARE_FECORE_CLASS();
};
