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
#include <FECore/FENLConstraint.h>
#include <FECore/FEEdge.h>
#include "febiomech_api.h"

class FEBIOMECH_API FESlideLine : public FEEdge
{
public:
	class FESlidingPoint
	{
	public:
		FESlidingPoint();

		void Serialize(DumpStream& ar);

		void Init();

	public:
		FELineElement*		pme = nullptr;
		double				gap = 0;
		vec3d				nu;
		double				r = 0;	//!< natural coordinates of primary surface projection on secondary surface element
		double				Ln = 0;	//!< contact pressure
		double				Lm = 0;	//!< Lagrange multipliers for contact pressure
	};

public:
	struct Projection
	{
		FELineElement* pe = nullptr;
		vec3d q;
		double r = 0;
	};

public:
	//! constructor
	FESlideLine(FEModel* pfem);

	//! Initializes data structures
	bool Init();

	FEMaterialPoint* CreateMaterialPoint() override;

	//! Serialize data to archive
	void Serialize(DumpStream& ar);

	Projection ClosestProjection(const vec3d& x);

public:
	vector<FESlidingPoint>		m_data;	//!< sliding contact surface data
};

class FEBIOMECH_API FESlideLineConstraint : public FENLConstraint
{
public:
	//! constructor
	FESlideLineConstraint(FEModel* pfem);

	//! Initializes sliding interface
	bool Init() override;

	//! interface activation
	void Activate() override;

	//! projects primary line nodes onto secondary line nodes
	void ProjectPoints(bool bupseg);

	//! calculate penalty value
	double Penalty() { return m_eps; }

	//! serialize data to archive
	void Serialize(DumpStream& ar) override;

	//! calculate contact pressures for file output
	void UpdateContactPressures();

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
	//! calculate the nodal force of a primary surface node
	void ContactNodalForce(int m, FELineElement& mel, vector<double>& fe);

	//! calculate the stiffness contribution of a single primary surface node
	void ContactNodalStiffness(int m, FELineElement& mel, matrix& ke);

private:
	void SerializePointers(DumpStream& ar);

public:
	int				m_laugon;
	int				m_naugmax;	//!< maximum nr of augmentations
	int				m_naugmin;	//!< minimum nr of augmentations
	double			m_gtol;		//!< gap tolerance
	double			m_atol;		//!< augmentation tolerance

	double			m_eps;		//!< penalty scale factor 

	int				m_nsegup;	//!< segment update parameter

private:
	bool	m_bfirst;	//!< flag to indicate the first time we enter Update
	double	m_normg0;	//!< initial gap norm

private:
	FESlideLine pl; // primary line
	FESlideLine sl; // secondary line

public:
	DECLARE_FECORE_CLASS();
};
