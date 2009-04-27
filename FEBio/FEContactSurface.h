#pragma once

#include "FESurface.h"
#include "vec2d.h"

//-----------------------------------------------------------------------------
//! This class describes a contact slave or master surface

//!	this class is used in contact analyses to describe a contacting
//! surface in a contact interface.

class FEContactSurface : public FESurface
{
public:
	//! constructor
	FEContactSurface(FEMesh* pm=0) : FESurface(pm) { m_NQ.Attach(this); }

	//! Initializes data structures
	void Init();

	//! shallow copy
	void ShallowCopy(FEContactSurface& s)
	{
		Lm  = s.Lm;
		gap = s.gap;
		pme.zero();
		Lt  = s.Lt;
	}

	//! Update the surface data
	void Update() {}

	//! Find element that contains the projection of x
	FEElement* FindMasterSegment(vec3d& x, vec3d& q, vec2d& r, bool& binit_nq, double tol);

	//! Calculate the total traction at a node
	vec3d traction(int inode);

public:
	vector<double>		gap;	//!< gap function at nodes
	vector<vec3d>		nu;		//!< master normal at slave node
	vector<FEElement*>	pme;	//!< master element a slave node penetrates
	vector<vec2d>		rs;		//!< natural coordinates of slave projection on master element
	vector<vec2d>		rsp;	//!< natural coordinates at previous time step
	vector<double>		Lm;		//!< Lagrange multipliers for contact pressure
	vector<mat2d>		M;		//!< surface metric tensor
	vector<vec2d>		Lt;		//!< Lagrange multipliers for friction
	vector<double>		off;	//!< gap offset (= shell thickness)
	vector<double>		eps;	//!< normal penalty factors

	FENNQuery		m_NQ;		//!< this structure is used in finding the master element that corresponds to a slave node
};
