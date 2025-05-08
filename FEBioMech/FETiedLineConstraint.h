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

class FETiedLine : public FEEdge
{
public:
	struct NodeData
	{
		FELineElement* me = nullptr;
		double r;
		vec3d vgap;
		vec3d Tc;
		vec3d Lm;
	};

	struct Projection
	{
		FELineElement* pe = nullptr;
		vec3d q;
		double r = 0;
	};

public:
	FETiedLine(FEModel* fem);

	FEMaterialPoint* CreateMaterialPoint() override;

	void Update();

	bool Create(FESegmentSet& eset) override;

	bool Init() override;

	Projection ClosestProjection(const vec3d& x);

public:
	std::vector<NodeData> m_data;
};


//! This class implements a tied-line.
class FETiedLineConstraint : public FENLConstraint
{
public:
	//! constructor
	FETiedLineConstraint(FEModel* pfem);

	//! destructor
	virtual ~FETiedLineConstraint(){}

	//! Initializes sliding interface
	bool Init() override;

	//! interface activation
	void Activate() override;

	//! serialize data to archive
	void Serialize(DumpStream& ar) override;

	//! build the matrix profile for use in the stiffness matrix
	void BuildMatrixProfile(FEGlobalMatrix& K) override;

	void ProjectLines(FETiedLine& pl, FETiedLine& sl);

public:
	//! calculate contact forces
	void LoadVector(FEGlobalVector& R, const FETimeInfo& tp) override;

	//! calculate contact stiffness
	void StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp) override;

	//! calculate Lagrangian augmentations
	bool Augment(int naug, const FETimeInfo& tp) override;

	//! update
	void Update() override;

public:
	int			m_laugon;	//!< enforcement method (0=penalty, 1=aug. Lag.)
	double		m_atol;		//!< augmentation tolerance
	double		m_eps;		//!< penalty scale factor
	double		m_stol;		//!< search tolerance
	int			m_naugmax;	//!< maximum nr of augmentations
	int			m_naugmin;	//!< minimum nr of augmentations
	bool		m_boffset;	//!< offset primary surface for shells
	double		m_Dmax;		//!< max distance for contact

private:
	FETiedLine pl; // primary line
	FETiedLine sl; // secondary line

	DECLARE_FECORE_CLASS();
};
