/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
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
#include "FEContactInterface.h"
#include "FEContactSurface.h"

//-----------------------------------------------------------------------------
//! This class describes a contact surface used for sticky contact.

//!	this class is used in contact analyses to describe a contacting
//! surface in a sticky contact interface.

class FEStickySurface : public FEContactSurface
{
public:
	class Data 
	{
	public:
		Data();

		void Serialize(DumpStream& ar);

	public:
		vec3d				gap;	//!< "gap" function
		vec2d				rs;		//!< natural coordinates of projection on secondary surface element
		vec3d				Lm;		//!< Lagrange multiplier
		vec3d				tn;		//!< traction vector
		FESurfaceElement*	pme;	//!< secondary surface element a node penetrates
	};

public:
	//! constructor
	FEStickySurface(FEModel* pfem) : FEContactSurface(pfem) {}

	//! Initializes data structures
	bool Init();

	//! data serialization
	void Serialize(DumpStream& ar);

public:
    void GetContactTraction(int nface, vec3d& pt);
	void GetNodalContactPressure(int nface, double* pn);
	void GetNodalContactTraction(int nface, vec3d* tn);

public:
	vector<Data>	m_data;	//!< node contact data
};

//-----------------------------------------------------------------------------
//! This class implements a sticky interface.
//! A sticky interface is like tied, but nodes are only tied when they come into
//! contact.

class FEStickyInterface : public FEContactInterface
{
public:
	//! constructor
	FEStickyInterface(FEModel* pfem);

	//! destructor
	virtual ~FEStickyInterface(){}

	//! Initializes sliding interface
	bool Init() override;

	//! interface activation
	void Activate() override;

	//! projects nodes onto secondary surface
	void ProjectSurface(FEStickySurface& ss, FEStickySurface& ms, bool bmove = false);

	//! serialize data to archive
	void Serialize(DumpStream& ar) override;

	//! return the primary and secondary surface
	FESurface* GetPrimarySurface() override { return &ss; }
	FESurface* GetSecondarySurface() override { return &ms; }

	//! return integration rule class
	bool UseNodalIntegration() override { return true; }

	//! build the matrix profile for use in the stiffness matrix
	void BuildMatrixProfile(FEGlobalMatrix& K) override;

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
	void SerializePointers(FEStickySurface& ss, FEStickySurface& ms, DumpStream& ar);

public:
	FEStickySurface	ss;	//!< primary surface
	FEStickySurface	ms;	//!< secondary surface

	double		m_atol;		//!< augmentation tolerance
	double		m_eps;		//!< penalty scale factor
	double		m_stol;		//!< search tolerance
	int			m_naugmax;	//!< maximum nr of augmentations
	int			m_naugmin;	//!< minimum nr of augmentations
	double		m_tmax;		//!< max traction
	double		m_snap;		//!< snap tolerance

	DECLARE_FECORE_CLASS();
};
