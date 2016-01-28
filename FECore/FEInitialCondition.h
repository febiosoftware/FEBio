#pragma once
#include "FEModelComponent.h"
#include <vector>
using namespace std;

//-----------------------------------------------------------------------------
//! Base class for defining initial conditions.
//! Initial conditions can be used to set the initial state of the model in an analysis. 
class FEInitialCondition : public FEModelComponent
{
public:
	FEInitialCondition(FEModel* pfem);
};

//-----------------------------------------------------------------------------
// Class representing an initial condition on a degree of freedom
class FEInitialBC : public FEInitialCondition
{
	struct ITEM
	{
		int		nid;	//!< node ID
		double	v;		//!< initial value
	};

public:
	FEInitialBC(FEModel* pfem);

	void SetDOF(int ndof) { m_dof = ndof; }

	void Serialize(DumpStream& ar);

	void Activate();

	void Add(int nid, double v) { ITEM it = {nid, v}; m_item.push_back(it); }

public:
	vector<ITEM>	m_item;		//!< node value pairs
	int				m_dof;		//!< degree of freedom
};

//-----------------------------------------------------------------------------
// Class for initializing degrees of freedom using a vec3d (useful for e.g. velocity)
class FEInitialBCVec3D : public FEInitialCondition
{
	struct ITEM
	{
		int		nid;	//!< node ID
		vec3d	v0;		//!< initial value
	};

public:
	FEInitialBCVec3D(FEModel* pfem) : FEInitialCondition(pfem) { m_dof[0] = m_dof[1] = m_dof[2] = -1; }

	void SetDOF(int d0, int d1, int d2) { m_dof[0] = d0; m_dof[1] = d1; m_dof[2] = d2; }

	void Serialize(DumpStream& ar);

	void Activate();

	void Add(int nid, vec3d v) { ITEM it = {nid, v}; m_item.push_back(it); }

public:
	vector<ITEM>	m_item;
	int				m_dof[3];
};
