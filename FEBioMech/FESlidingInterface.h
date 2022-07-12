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
#include "FEContactSurface.h"
#include "FEContactInterface.h"
#include <FECore/FEClosestPointProjection.h>
#include <FECore/vector.h>

//-----------------------------------------------------------------------------
class FEBIOMECH_API FESlidingSurface : public FEContactSurface
{
	class FESlidingPoint : public FEContactMaterialPoint 
	{
	public:
		FESlidingPoint();

		void Serialize(DumpStream& ar) override;

		void Init() override;

	public:
		vec3d				m_nu;	  //!< secondary surface normal at primary surface node
		vec2d				m_rs;	  //!< natural coordinates of primary surface projection on secondary surface element
		vec2d				m_rsp;  //!< natural coordinates at previous time step
		double				m_Lm;	  //!< Lagrange multipliers for contact pressure
		mat2d				m_M;	  //!< surface metric tensor
		vec2d				m_Lt;	  //!< Lagrange multipliers for friction
		double				m_off;  //!< gap offset (= shell thickness)
		double				m_eps;  //!< normal penalty factors
	};

public:
	//! constructor
	FESlidingSurface(FEModel* pfem);

	//! Initializes data structures
	bool Init();

	//! Calculate the total traction at a node
	vec3d traction(int inode);

	//! evaluate net contact force
	vec3d GetContactForce();
	
	//! evaluate net contact area
	double GetContactArea();
    
	//! Serialize data to archive
	void Serialize(DumpStream& ar);

public:
    void GetContactTraction(int nface, vec3d& pt);
	void GetNodalContactPressure(int nface, double* pg);
	void GetNodalContactTraction(int nface, vec3d* pt);

public:
	vector<FESlidingPoint>		m_data;	//!< sliding contact surface data
};

//-----------------------------------------------------------------------------
//! This class implements a sliding interface

//! The FESlidingInterface class defines an interface between two surfaces.

class FEBIOMECH_API FESlidingInterface : public FEContactInterface
{
public:
	//! constructor
	FESlidingInterface(FEModel* pfem);

	//! destructor
	virtual ~FESlidingInterface(){}

	//! Initializes sliding interface
	bool Init() override;

	//! interface activation
	void Activate() override;

	//! projects primary surface nodes onto secondary surface nodes
	void ProjectSurface(FESlidingSurface& ss, FESlidingSurface& ms, bool bupseg, bool bmove = false);

	//! calculate penalty value
	double Penalty() { return m_eps; }

	//! serialize data to archive
	void Serialize(DumpStream& ar) override;

	//! calculate contact pressures for file output
	void UpdateContactPressures();

	//! return the primary and secondary surface
	FESurface* GetPrimarySurface() override { return &m_ss; }
	FESurface* GetSecondarySurface() override { return &m_ms; }

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

protected:
	//! calculate auto penalty factor
	void CalcAutoPenalty(FESlidingSurface& s);

	//! calculate the nodal force of a primary surface node
	void ContactNodalForce(int m, FESlidingSurface& ss, FESurfaceElement& mel, vector<double>& fe);

	//! calculate the stiffness contribution of a single primary surface node
	void ContactNodalStiffness(int m, FESlidingSurface& ss, FESurfaceElement& mel, matrix& ke);

	//! map the frictional data from the old element to the new element
	void MapFrictionData(int inode, FESlidingSurface& ss, FESlidingSurface& ms, FESurfaceElement& sn, FESurfaceElement& so, vec3d& q);

private:
	void SerializePointers(FESlidingSurface& ss, FESlidingSurface& ms, DumpStream& ar);

public:
	FESlidingSurface	m_ss;	//!< primary surface
	FESlidingSurface	m_ms;	//!< secondary surface

	bool			m_btwo_pass;	//!< two pass algorithm flag

	double			m_sradius;			//!< search radius for self contact

	int				m_naugmax;	//!< maximum nr of augmentations
	int				m_naugmin;	//!< minimum nr of augmentations
	double			m_gtol;		//!< gap tolerance
	double			m_atol;		//!< augmentation tolerance
	double			m_ktmult;	//!< multiplier for tangential stiffness
	double			m_knmult;	//!< multiplier for normal stiffness

	double			m_stol;		//!< search tolerance

	bool			m_bautopen;	//!< auto penalty calculation
	double			m_eps;		//!< penalty scale factor 
	bool			m_bupdtpen;	//!< update penalty 

	bool			m_breloc;	//!< initial node relocation

	double			m_mu;		//!< friction coefficient
	double			m_epsf;		//!< penalty scale factor for friction

	int				m_nsegup;	//!< segment update parameter

private:
	bool	m_bfirst;	//!< flag to indicate the first time we enter Update
	double	m_normg0;	//!< initial gap norm

public:
	DECLARE_FECORE_CLASS();
};
