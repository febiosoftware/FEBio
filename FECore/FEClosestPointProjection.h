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
	FESurfaceElement* Project(vec3d& x, vec3d& q, vec2d& r, double tol);

protected:
	FESurface&		m_surf;		//!< reference to surface
	FENNQuery		m_SNQ;		//!< used to find the nearest neighbour
	FENodeElemTree	m_NET;		//!< node-element tree
};
