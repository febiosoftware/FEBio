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



#include "stdafx.h"
#include "FEClosestPointProjection.h"
#include "FEElemElemList.h"
#include "FEMesh.h"

//-----------------------------------------------------------------------------
// constructor
FEClosestPointProjection::FEClosestPointProjection(FESurface& s) : m_surf(s)
{
	// set default options
	m_tol = 0.01;
	m_rad = 0.0;	// 0 means don't use search radius
	m_bspecial = false;
	m_projectBoundary = false;

	// calculate node-element list
	m_NEL.Create(m_surf);
	m_EEL.Create(&m_surf);
}

//-----------------------------------------------------------------------------
//! Initialization of data structures
bool FEClosestPointProjection::Init()
{
	// initialize the nearest neighbor search
	m_SNQ.Attach(&m_surf);
	m_SNQ.Init();

	return true;
}

//-----------------------------------------------------------------------------
// helper function for projecting a point onto an edge
bool Project2Edge(const vec3d& p0, const vec3d& p1, const vec3d& x, vec3d& q)
{
	vec3d e = p1 - p0;
	double D = e*e;
	if (D == 0.0) return false;

	double l = (e*(x - p0))/D;

	q = p0 + e*l;

	return (l>=0.0)&&(l<=1.0);
}


//-----------------------------------------------------------------------------
class FEPatch
{
public:
	FEPatch(FESurface* surf, FEElement** ppe, int n)
	{
		m_surf = surf;
		m_ppe = ppe;
		m_nval = n;
	}

	bool HasNode(int n)
	{
		for (int i = 0; i < m_nval; ++i)
		{
			// get the element
			FESurfaceElement& el = static_cast<FESurfaceElement&> (*m_ppe[i]);
			if (el.HasNode(n))
			{
				return true;
			}
		}
		return false;
	}

	bool Contains(FESurfaceElement& el)
	{
		for (int j = 0; j < m_nval; ++j)
		{
			FESurfaceElement& ei = *Element(j);
			if ((el.Type() == ei.Type()) && (el.Nodes() == ei.Nodes()))
			{
				if (ei.HasNodes(&(el.m_node)[0], el.Nodes()) != 0)
				{
					return true;
				}
			}
		}
		return false;
	}

	FESurface* GetSurface() { return m_surf; }

	int Size() const { return m_nval; }
	FESurfaceElement* Element(int n) { return dynamic_cast<FESurfaceElement*>(m_ppe[n]); }

private:
	FESurface*	m_surf;
	FEElement**	m_ppe;
	int			m_nval;
};

// project a point to a patch
FESurfaceElement* projectToPatch(FEPatch& patch, const vec3d& x, vec3d& q, double& r, double& s, double tol)
{
	FESurfaceElement* pemin = nullptr;
	int nval = patch.Size();
	FESurface* surf = patch.GetSurface();
	double d2min = 0;
	for (int j = 0; j < nval; ++j)
	{
		// get the element
		FESurfaceElement& el = *patch.Element(j);
		int N = el.Nodes();

		// project the node on the element
		double p[2] = { 0, 0 };
		vec3d qj = surf->ProjectToSurface(el, x, p[0], p[1]);
		double d2 = (qj - x).norm2();
		if (surf->IsInsideElement(el, p[0], p[1], tol))
		{
			if ((pemin == nullptr) || (d2 < d2min))
			{
				d2min = d2;
				q = qj;
				r = p[0];
				s = p[1];
				pemin = &el;
			}
		}
	}

	return pemin;
}

//-----------------------------------------------------------------------------
//! Project the node on the surface. This function returns a pointer to the
//! surface element if the projection is succesful, otherwise null. The natural
//! coordinates of the projection is return in r and the spatial coordinates in q.
//! 
FESurfaceElement* FEClosestPointProjection::Project(const vec3d& x, vec3d& q, vec2d& r)
{
	// get the mesh
	FEMesh& mesh = *m_surf.GetMesh();

	// let's find the closest node
	int mn = m_SNQ.Find(x);
	if (mn < 0) return nullptr;

	// make sure it is within the search radius
	vec3d rm = m_surf.Node(mn).m_rt;
	double d2 = (x - rm)*(x - rm);
	if ((m_rad > 0) && (d2 > m_rad))
	{
		return nullptr;
	}

	// now that we found the closest node, lets see if we can find 
	// the best element
	FEPatch patch(&m_surf, m_NEL.ElementList(mn), m_NEL.Valence(mn));
	FESurfaceElement* pe = projectToPatch(patch, x, q, r[0], r[1], m_tol);
	if (pe) return pe;

	// If we get here, we did not find a facet.
	// There are a couple of reasons why the search could fail:
	// -1. the point cannot be projected onto the surface. For contact this implies the node is not in contact.
	// -2. the projection falls outside the set of elements surrounding the closest point.
	// -3. the projection falls on an edge of two faces whos normals are pointing away.
	// -4. the closest node is in fact the closest point and no closer projection on face or edge can be found
	//
	if (m_bspecial)
	{
		return ProjectSpecial(mn, x, q, r);
	}

	return nullptr;
}

//-----------------------------------------------------------------------------
//! This is slightly modified version of the ClosestPointProjection where the
//! node to be projected is part of this surface. This is used in some contact algorithms
//! that support self-contact. The algorithm is modified in two ways. First, a search_radius
//! limits the max distance for consideration of closest point. Second, the search point cannot
//! be part of the star of the closest point.
FESurfaceElement* FEClosestPointProjection::Project(int nodeIndex, vec3d& q, vec2d& r)
{
	// get the mesh
	FEMesh& mesh = *m_surf.GetMesh();

	// get the node's position
	vec3d x = mesh.Node(nodeIndex).m_rt;

	// Find the closest surface node to x that:
	// 1. is within the search radius
	// 2. its star does not contain n
	int mn = -1;	// local index of closest node
	double d2min;	// min squared distance
	int N = m_surf.Nodes();
	double R2 = m_rad * m_rad;
	for (int i = 0; i < N; ++i)
	{
		if (m_surf.NodeIndex(i) != nodeIndex)
		{
			vec3d r = m_surf.Node(i).m_rt;
			double d2 = (r - x)*(r - x);

			// make sure the node lies within the search radius
			if ((m_rad == 0) || (d2 <= R2))
			{
				bool bok = true;

				// The node cannot be part of the star of the closest point
				FEPatch patch(&m_surf, m_NEL.ElementList(i), m_NEL.Valence(i));
				if (patch.HasNode(nodeIndex))
				{
					bok = false;
				}

				// make sure the point is closer than the last one
				if (bok && ((mn == -1) || (d2 < d2min)))
				{
					d2min = d2;
					mn = i;
					q = r;
				}
			}
		}
	}
	if (mn == -1) return nullptr;

	// now that we found the closest node, lets see if we can find 
	// the best element
	FEPatch patch(&m_surf, m_NEL.ElementList(mn), m_NEL.Valence(mn));
	FESurfaceElement* pemin = projectToPatch(patch, x, q, r[0], r[1], m_tol);
	if (pemin) return pemin;

	// If we get here, we did not find a facet.
	// There are a couple of reasons why the search has failed:
	// -1. the point cannot be projected onto the surface. For contact this implies the node is not in contact.
	// -2. the projection falls outside the set of elements surrounding the closest point.
	// -3. the projection falls on an edge of two faces whos normals are pointing away.
	// -4. the closest node is in fact the closest point and no closer projection on face or edge can be found
	//
	if (m_bspecial)
	{
		return ProjectSpecial(mn, x, q, r);
	}

	return nullptr;
}

//! Project a point of a surface element onto a surface
FESurfaceElement* FEClosestPointProjection::Project(FESurfaceElement* pse, int intgrPoint, vec3d& q, vec2d& r)
{
	assert(pse);
	if (pse == nullptr) return nullptr;

	// get the mesh
	FEMesh& mesh = *m_surf.GetMesh();

	// get the node's position
	int nn = pse->Nodes();
	vec3d re[FEElement::MAX_NODES];
	for (int l = 0; l < nn; ++l) re[l] = mesh.Node(pse->m_node[l]).m_rt;
	vec3d x = pse->eval(re, intgrPoint);

	// if the element belongs to this surface
	// we want to prevent this element and its immediate neigbors
	// from being the closest.
	bool check_self_projection = false;
	if (ContainsElement(pse))
	{
		check_self_projection = true;
	}

	// find the closest point
	int mn = -1;
	double d2min;
	double R2 = m_rad * m_rad;
	int N = m_surf.Nodes();
	for (int i = 0; i < N; ++i)
	{
		vec3d ri = m_surf.Node(i).m_rt;
		double d2 = (ri - x)*(ri - x);

		// make sure the node lies within the search radius
		if ((R2 == 0) || (d2 <= R2))
		{
			bool bok = true;
			if (check_self_projection)
			{
				// The pse element cannot be part of the star of the closest point
				FEPatch patch(&m_surf, m_NEL.ElementList(i), m_NEL.Valence(i));
				if (patch.Contains(*pse))
				{
					bok = false;
				}
			}

			if (bok && ((mn == -1) || (d2 < d2min)))
			{
				d2min = d2;
				mn = i;
				q = ri;
			}
		}
	}
	if (mn == -1) return nullptr;

	// mn is a local index, so get the global node number too
	int m = m_surf.NodeIndex(mn);

	// get the nodal position
	vec3d r0 = mesh.Node(m).m_rt;

	// check the distance
	double D = (x - r0).norm();
	if ((m_rad > 0) && (D > m_rad)) return 0;

	// now that we found the closest node, lets see if we can find 
	// the best element
	FEPatch patch(&m_surf, m_NEL.ElementList(mn), m_NEL.Valence(mn));
	FESurfaceElement* pe = projectToPatch(patch, x, q, r[0], r[1], m_tol);
	if (pe) return pe;

	// If we get here, we did not find a facet.
	// handle special cases
	if (m_bspecial)
	{
		return ProjectSpecial(mn, x, q, r);
	}

	return nullptr;
}

bool FEClosestPointProjection::ContainsElement(FESurfaceElement* el)
{
	if (el == nullptr) return false;

	int n = el->m_lnode[0];
	if ((n < 0) || (n >= m_surf.Nodes())) return false;

	int nval = m_NEL.Valence(n);
	FEElement** pe = m_NEL.ElementList(n);
	for (int i = 0; i < nval; ++i)
	{
		FESurfaceElement& ei = dynamic_cast<FESurfaceElement&>(*(pe[i]));

		int m = el->Nodes();
		if ((el->Type() == ei.Type()) && (ei.Nodes() == m))
		{
			if (el->HasNodes(&ei.m_node[0], m) != 0)
			{
				return true;
			}
		}
	}

	return false;
}

//-------------------------------------------------------------------------------------------
// This function handles special cases for closest-point searches
// -1. the point cannot be projected onto the surface. For contact this implies the node is not in contact.
// -2. the projection falls outside the set of elements surrounding the closest point.
// -3. the projection falls on an edge of two faces whos normals are pointing away.
// -4. the closest node is in fact the closest point and no closer projection on face or edge can be found
//
FESurfaceElement* FEClosestPointProjection::ProjectSpecial(int closestPoint, const vec3d& x, vec3d& q, vec2d& r)
{
	FEMesh& mesh = *m_surf.GetMesh();
	int nval = m_NEL.Valence(closestPoint);
	FEElement** pe = m_NEL.ElementList(closestPoint);
	int* eil = m_NEL.ElementIndexList(closestPoint);
	int m = m_surf.NodeIndex(closestPoint);

	// get position of closest point
	vec3d rm = m_surf.Node(closestPoint).m_rt;

	FESurfaceElement* pemin = nullptr;
	double D2min = (x - rm).norm2();

	// First, let's redo the search but with a larger search radius
	// NOTE: This is highly inefficient since multiple nodes and faces are visited multiple times
	for (int i = 0; i < nval; ++i)
	{
		// get the element
		FESurfaceElement& el = static_cast<FESurfaceElement&> (*pe[i]);
		int N = el.Nodes();
		for (int j = 0; j < N; ++j)
		{
			int mj = el.m_lnode[j];
			int nj = m_NEL.Valence(mj);
			FEElement** pej = m_NEL.ElementList(mj);
			for (int k = 0; k < nj; ++k)
			{
				FESurfaceElement& ek = static_cast<FESurfaceElement&> (*pej[k]);

				// project the node on the element
				double p[2] = { 0, 0 };
				vec3d qk = m_surf.ProjectToSurface(ek, x, p[0], p[1]);
				if (m_surf.IsInsideElement(ek, p[0], p[1], m_tol))
				{
					double D2 = (x - qk).norm2();
					if (D2 < D2min)
					{
						pemin = &ek;
						q = qk;
						r = vec2d(p[0], p[1]);
						D2min = D2;
					}
				}
			}
		}
	}
	if (pemin) return pemin;

	// This tries to handle some of the special cases. 
	// First we try to find an edge this point projects on
	for (int j = 0; j < nval; ++j)
	{
		// get an element
		FESurfaceElement& el = static_cast<FESurfaceElement&> (*pe[j]);

		// loop over facet edges
		int N = el.facet_edges();
		for (int k = 0; k < N; ++k)
		{
			// if projection on boundary edges are not allowed
			// make sure the element has a neighbor
			if (m_projectBoundary || m_EEL.Neighbor(eil[j], k))
			{
				// get the two edge node indices
				int en[3];
				el.facet_edge(k, en);
				int nk1 = en[0];
				int nk2 = en[1];

				// make sure one of them is our closest point
				if ((nk1 == closestPoint) || (nk2 == closestPoint))
				{
					// try to project it on the edge
					vec3d p0 = m_surf.Node(nk1).m_rt;
					vec3d p1 = m_surf.Node(nk2).m_rt;
					vec3d qk;
					if (Project2Edge(p0, p1, x, qk))
					{
						double p[2] = { 0, 0 };
						vec3d qp = m_surf.ProjectToSurface(el, qk, p[0], p[1]);
						if (m_surf.IsInsideElement(el, p[0], p[1], m_tol))
						{
							// see if this is a closer projection
							double D2 = (x - qp).norm2();
							if (D2 < D2min)
							{
								pemin = &el;
								q = qp;
								D2min = D2;
								r = vec2d(p[0], p[1]);
							}
						}
					}
				}
			}
		}
	}
	if (pemin) return pemin;

	// if we get here then no edge was found.
	// This can imply that the projection is in fact the closest node

	// first we want to make sure that the node does not lie on the boundary if
	// boundary projects are not allowed
	if (m_projectBoundary == false)
	{
		for (int j = 0; j < nval; ++j)
		{
			// get an element
			FESurfaceElement& el = static_cast<FESurfaceElement&> (*pe[j]);

			// make sure that the node does not lie on a boundary edge
			int N = el.facet_edges();
			for (int k = 0; k < N; ++k)
			{
				if (m_EEL.Neighbor(eil[j], k) == 0)
				{
					int n0 = el.m_node[k];
					int n1 = el.m_node[(k + 1) % N];
					if ((n0 == m) || (n1 == m)) return 0;
				}
			}
		}
	}

	// Any element can now be used but just to be safe, we loop over all
	// and make sure the element actually contains the closest point.
	for (int j = 0; j < nval; ++j)
	{
		// get an element
		FESurfaceElement& el = static_cast<FESurfaceElement&> (*pe[j]);

		// make sure that the node does not lie on a boundary edge
		int N = el.facet_edges();
		if (m_projectBoundary == false)
		{
			for (int k = 0; k < N; ++k)
			{
				if (m_EEL.Neighbor(eil[j], k) == 0)
				{
					int n0 = el.m_node[k];
					int n1 = el.m_node[(k + 1) % N];
					if ((n0 == m) || (n1 == m)) return 0;
				}
			}
		}

		// make sure one of them is our closest point
		if (el.HasNode(m))
		{
			q = rm;
			return &el;
		}
	}

	return 0;
}
