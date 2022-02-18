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
#include "FECore/FEModel.h"
#include <FECore/tens4d.h>

//-----------------------------------------------------------------------------
// Class describing the RVE model.
// This is used by the homogenization code.
class FERVEModel : public FEModel
{
public:
	enum RVE_TYPE
	{
		DISPLACEMENT,		// prescribed displacement
		PERIODIC_LC,		// periodic, linear constraints
		PERIODIC_AL			// periodic, augmented Lagrangian (NOTE: obsolete, should probably delete)
	};

public:
	FERVEModel();
	~FERVEModel();

	//! initialization
	bool Init() override;

	//! one time initialization (only called on the master RVE)
	bool InitRVE(int rveType, const char* szbc);

	//! Return the initial volume (calculated in Init)
	double InitialVolume() const { return m_V0; }

	//! return current volume (calculated each time)
	double CurrentVolume();

	// scale the geometry
	void ScaleGeometry(double scale);

	//! see if node is boundary node
	bool IsBoundaryNode(int i) const { return (m_BN[i]==1); }

	//! Update the RVE (before it is solved)
	void Update(const mat3d& F);

	// copy from the master RVE
	void CopyFrom(FERVEModel& rve);

	// set the parent FEModel
	void SetParentModel(FEModel* fem);

	//! Calculate the stress average
	mat3ds StressAverage(mat3d& F, FEMaterialPoint& mp);
	mat3ds StressAverage(FEMaterialPoint& mp);

	//! Calculate the stiffness average
	tens4ds StiffnessAverage(FEMaterialPoint &mp);

protected:
	//! Calculate the initial volume
	void EvalInitialVolume();

	//! find the list of boundary nodes
	void FindBoundaryNodes(vector<int>& BN);

	//! Center the RVE
	void CenterRVE();

	bool PrepDisplacementBC(FENodeSet* set);
	bool PrepPeriodicBC(const char* szbc);
	bool PrepPeriodicLC();

private:
	// Hide the solve function so we can't call it. 
	// (We use the RCI solution method)
	bool Solve() override;

protected:
	FEModel*		m_parentfem;		//!< parent FEModel
	double			m_V0;				//!< initial volume
	int				m_bctype;			//!< RVE type
	FEBoundingBox	m_bb;				//!< bounding box of mesh
	vector<int>		m_BN;				//!< boundary node flags
};
