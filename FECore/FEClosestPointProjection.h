#pragma once
#include "FESurface.h"
#include "FENNQuery.h"

//-----------------------------------------------------------------------------
// This class can be used to find the closest point projection of a point
// onto a surface.
class FEClosestPointProjection
{
public:
	//! constructor
	FEClosestPointProjection(FESurface& s);

	//! Initialization
	bool Init();

	//! Project point onto surface
	FESurfaceElement* Project(vec3d& x, vec3d& q, vec2d& r);

public:
	//! Set the projection tolerance
	void SetTolerance(double t) { m_tol = t; }

	//! get the projection tolerance
	double GetTolerance() { return m_tol; }

protected:
	double	m_tol;	//!< projection tolerance

protected:
	FESurface&		m_surf;		//!< reference to surface
	FENNQuery		m_SNQ;		//!< used to find the nearest neighbour
	FENodeElemTree	m_NET;		//!< node-element tree
};
