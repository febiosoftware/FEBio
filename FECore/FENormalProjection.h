#pragma once
#include "FESurface.h"
#include "FEOctree.h"

//-----------------------------------------------------------------------------
//! This class calculates the normal projection on to a surface.
//! This is used by some contact algorithms.
class FENormalProjection
{
public:
	//! constructor
	FENormalProjection(FESurface& s);

	// initialization
	void Init();

	void SetTolerance(double tol) { m_tol = tol; }
	void SetSearchRadius(double srad) { m_rad = srad; }

public:
	//! find the intersection of a ray with the surface
	FESurfaceElement* Project(vec3d r, vec3d n, double rs[2]);
	FESurfaceElement* Project2(vec3d r, vec3d n, double rs[2]);
	FESurfaceElement* Project3(const vec3d& r, const vec3d& n, double rs[2], int* pei = 0);

	vec3d Project(const vec3d& r, const vec3d& N);

private:
	double	m_tol;	//!< projection tolerance
	double	m_rad;	//!< search radius

private:
	FESurface&	m_surf;	//!< the target surface
	FEOctree	m_OT;	//!< used to optimize ray-surface intersections
};
