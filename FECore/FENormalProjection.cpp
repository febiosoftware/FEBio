#include "stdafx.h"
#include "FENormalProjection.h"
#include "FEMesh.h"

//-----------------------------------------------------------------------------
FENormalProjection::FENormalProjection(FESurface& s) : m_surf(s)
{
	m_tol = 0.0;
	m_rad = 0.0;
}

//-----------------------------------------------------------------------------
void FENormalProjection::Init()
{
	m_OT.Attach(&m_surf);
	m_OT.Init(m_tol);
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

//-----------------------------------------------------------------------------
//! This function finds the element which is intersected by the ray (r,n).
//! It returns a pointer to the element, as well as the isoparametric coordinates
//! of the intersection point.
//!
FESurfaceElement* FENormalProjection::Project3(const vec3d& r, const vec3d& n, double rs[2], int* pei)
{
	// let's find all the candidate surface elements
	set<int>selist;
	m_OT.FindCandidateSurfaceElements(r, n, selist);

	double g, gmax = -1e99, r2[2] = {rs[0], rs[1]};
	int imin = -1;
	FESurfaceElement* pme = 0;

	// loop over all surface element
	set<int>::iterator it;
	for (it = selist.begin(); it != selist.end(); ++it)
	{
		FESurfaceElement& el = m_surf.Element(*it);

		// see if the ray intersects this element
		if (m_surf.Intersect(el, r, n, r2, g, m_tol))
		{
			// see if this is the best intersection found so far
			// TODO: should I put a limit on how small g can
			//       be to be considered a valid intersection?
//			if (g < gmin)
			if (g > gmax)
			{
				// keep results
				pme = &el;
//				gmin = g;
				gmax = g;
				imin = *it;
				rs[0] = r2[0];
				rs[1] = r2[1];
			}
		}	
	}

	if (pei) *pei = imin;

	// return the intersected element (or zero if none)
	return pme;
}

//-----------------------------------------------------------------------------
vec3d FENormalProjection::Project(const vec3d& x, const vec3d& N)
{
	double rs[2];
	FESurfaceElement* pe = FENormalProjection::Project3(x, N, rs);
	if (pe)
	{
		FEMesh& mesh = *m_surf.GetMesh();
		vec3d r[FEElement::MAX_NODES];
		for (int i=0; i<pe->Nodes(); ++i) r[i] = mesh.Node(pe->m_node[i]).m_rt;
		vec3d q = pe->eval(r, rs[0], rs[1]);
		return q;
	}
	else return x;
}
