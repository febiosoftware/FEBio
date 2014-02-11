#include "stdafx.h"
#include "FENormalProjection.h"

//-----------------------------------------------------------------------------
FENormalProjection::FENormalProjection(FESurface& s) : m_surf(s)
{
}

//-----------------------------------------------------------------------------
void FENormalProjection::Init()
{
	m_OT.Attach(&m_surf);
	m_OT.Init();
}

//-----------------------------------------------------------------------------
//! This function finds the element which is intersected by the ray (r,n).
//! It returns a pointer to the element, as well as the isoparametric coordinates
//! of the intersection point. It searches for the closest patch based on
//! algebraic value of the gap function
//!
FESurfaceElement* FENormalProjection::Project(vec3d r, vec3d n, double rs[2])
{
	// let's find all the candidate surface elements
	set<int>selist;
	m_OT.FindCandidateSurfaceElements(r, n, selist);
	
	// now that we found candidate surface elements, lets see if we can find 
	// those that intersect the ray, then pick the closest intersection
	set<int>::iterator it;
	bool found = false;
	double rsl[2], gl, g;
	FESurfaceElement* pei = 0;
	for (it=selist.begin(); it!=selist.end(); ++it) {
		// get the surface element
		int j = *it;
		// project the node on the element
		FESurfaceElement* pe = &m_surf.Element(j);
		if (m_surf.Intersect(*pe, r, n, rsl, gl, m_tol)) {
			if ((!found) && (gl > -m_rad)) {
				found = true;
				g = gl;
				rs[0] = rsl[0];
				rs[1] = rsl[1];
				pei = pe;
			} else if ((gl < g) && (gl > -m_rad)) {
				g = gl;
				rs[0] = rsl[0];
				rs[1] = rsl[1];
				pei = pe;
			}
		}
	}
	if (found) return pei;
	
	// we did not find a master surface
	return 0;
}

//-----------------------------------------------------------------------------
//! This function finds the element which is intersected by the ray (r,n).
//! It returns a pointer to the element, as well as the isoparametric coordinates
//! of the intersection point.  It searches for the closest patch based on
//! the absolute value of the gap function
//!
FESurfaceElement* FENormalProjection::Project2(vec3d r, vec3d n, double rs[2])
{
	// let's find all the candidate surface elements
	set<int>selist;
	m_OT.FindCandidateSurfaceElements(r, n, selist);
	
	// now that we found candidate surface elements, lets see if we can find 
	// those that intersect the ray, then pick the closest intersection
	set<int>::iterator it;
	bool found = false;
	double rsl[2], gl, g;
	FESurfaceElement* pei = 0;
	for (it=selist.begin(); it!=selist.end(); ++it) {
		// get the surface element
		int j = *it;
		FESurfaceElement* pe = &m_surf.Element(j);
		// project the node on the element
		if (m_surf.Intersect(*pe, r, n, rsl, gl, m_tol)) {
			if ((!found) && (fabs(gl) < m_rad)) {
				found = true;
				g = gl;
				rs[0] = rsl[0];
				rs[1] = rsl[1];
				pei = pe;
			} else if ((fabs(gl) < fabs(g)) && (fabs(gl) < m_rad)) {
				g = gl;
				rs[0] = rsl[0];
				rs[1] = rsl[1];
				pei = pe;
			}
		}
	}
	if (found) return pei;
	
	// we did not find a master surface
	return 0;
}

