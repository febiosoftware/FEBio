/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2020 University of Utah, The Trustees of Columbia University in 
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
	m_handleQuads = false;

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
//! Project the node on the surface. This function returns a pointer to the
//! surface element if the projection is succesful, otherwise null. The natural
//! coordinates of the projection is return in r and the spatial coordinates in q.
//! 
FESurfaceElement* FEClosestPointProjection::Project(vec3d& x, vec3d& q, vec2d& r)
{
	// get the mesh
	FEMesh& mesh = *m_surf.GetMesh();

	// let's find the closest node
	int mn = m_SNQ.Find(x);

	// mn is a local index, so get the global node number too
	int m = m_surf.NodeIndex(mn);

	// get the nodal position
	vec3d rm = mesh.Node(m).m_rt;

	// now that we found the closest node, lets see if we can find 
	// the best element
	int nval = m_NEL.Valence(mn);
	FEElement** pe = m_NEL.ElementList(mn);
	int* eil = m_NEL.ElementIndexList(mn);
	for (int j=0; j<nval; ++j)
	{
		// get the element
		FESurfaceElement& el = static_cast<FESurfaceElement&> (*pe[j]);
		int N = el.Nodes();

		// project the node on the element
		r[0] = 0;
		r[1] = 0;
		q = m_surf.ProjectToSurface(el, x, r[0], r[1]);
		if (m_surf.IsInsideElement(el, r[0], r[1], m_tol)) return &el;
	}

	// If we get here, we did not find a facet.
	// There are a couple of reasons why the search has failed:
	// -1. the point cannot be projected onto the surface. For contact this implies the node is not in contact.
	// -2. the projection falls outside the set of elements surrounding the closest point.
	// -3. the projection falls on an edge of two faces whos normals are pointing away.
	// -4. the closest node is in fact the closest point and no closer projection on face or edge can be found
	//
	if (m_bspecial)
	{
		// First, let's redo the search but with a larger search radius
		// NOTE: This is highly inefficient since multiple nodes and faces are visited multiple times
		for (int i=0; i<nval; ++i)
		{
			// get the element
			FESurfaceElement& el = static_cast<FESurfaceElement&> (*pe[i]);
			int N = el.Nodes();
			for (int j=0; j<N; ++j)
			{
				int mj = el.m_lnode[j];
				int nj = m_NEL.Valence(mj);
				FEElement** pej = m_NEL.ElementList(mj);
				for (int k=0; k<nj; ++k)
				{
					FESurfaceElement& ek = static_cast<FESurfaceElement&> (*pej[k]);

					// project the node on the element
					r[0] = 0;
					r[1] = 0;
					q = m_surf.ProjectToSurface(ek, x, r[0], r[1]);
					if (m_surf.IsInsideElement(ek, r[0], r[1], m_tol)) return &ek;
				}
			}
		}

		// This tries to handle some of the special cases. 
		// For now, we only consider triangular facets. 
		// First we try to find an edge this point projects on
		for (int j=0; j<nval; ++j)
		{
			// get an element
			FESurfaceElement& el = static_cast<FESurfaceElement&> (*pe[j]);

			// only consider triangles (for now)
			int N = el.Nodes();

			// TODO: The case N==4 is causing convergence issues for some test suite problems.
			//       I'm commenting it out until I get a chance to look further into this
			if ((N == 3) || (m_handleQuads && (N == 4)))
			{
				int kmin = -1;
				// we got to find a point closer than the closest point
				double Dmin = (x - rm).norm();
				vec3d qmin;
				for (int k=0; k<N; ++k)
				{
					// if projection on boundary edges are not allowed
					// make sure the element has a neighbor
					if (m_projectBoundary || m_EEL.Neighbor(eil[j], k))
					{
						// get the two edge node indices
						int nk1 = el.m_node[k];
						int nk2 = el.m_node[(k+1)%N];

						// make sure one of them is our closest point
						if ((nk1 == m) || (nk2 == m))
						{
							// try to project it on the edge
							vec3d p0 = mesh.Node(nk1).m_rt;
							vec3d p1 = mesh.Node(nk2).m_rt;
							if (Project2Edge(p0, p1, x, q))
							{
								// see if this is a closer projection
								double D = (x-q).norm();
								if (D <= Dmin)
								{
									kmin = k;
									qmin = q;
									Dmin = D;
								}
							}
						}
					}
				}

				if (kmin != -1)
				{
					// okay, we got an edge.
					// Now find the iso-parametric coordinates
					q = m_surf.ProjectToSurface(el, qmin, r[0], r[1]);
					if (m_surf.IsInsideElement(el, r[0], r[1], m_tol)) return &el;
				}
			}
		}

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

				// only consider triangles (for now)
				int N = el.Nodes();
				if ((N == 3) || (m_handleQuads && (N==4)))
				{
					// make sure that the node does not lie on a boundary edge
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
		}

		// Any element can now be used but just to be safe, we loop over all
		// and make sure the element actually contains the closest point.
		for (int j=0; j<nval; ++j)
		{
			// get an element
			FESurfaceElement& el = static_cast<FESurfaceElement&> (*pe[j]);

			// only consider triangles (for now)
			int N = el.Nodes();
			if ((N == 3) || (m_handleQuads && (N==4)))
			{
				// make sure that the node does not lie on a boundary edge
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

				int* en = &(el.m_node[0]);

				// make sure one of them is our closest point
				if (N == 3)
				{
					q = rm;
					if (en[0] == m) { r[0] = 0; r[1] = 0; }
					if (en[1] == m) { r[0] = 1; r[1] = 0; }
					if (en[2] == m) { r[0] = 0; r[1] = 1; }
					return &el;
				}
				else if (N == 4)
				{
					q = rm;
					if (en[0] == m) { r[0] = -1; r[1] = -1; }
					if (en[1] == m) { r[0] =  1; r[1] = -1; }
					if (en[2] == m) { r[0] =  1; r[1] =  1; }
					if (en[3] == m) { r[0] = -1; r[1] =  1; }
					return &el;
				}
			}
		}
	}

	return 0;
}

//-----------------------------------------------------------------------------
//! This is slightly modified version of the ClosestPointProjection where the
//! node to be projected is part of this surface. This is used in some contact algorithms
//! that support self-contact. The algorithm is modified in two ways. First, a search_radius
//! limits the max distance for consideration of closest point. Second, the search point cannot
//! be part of the net of the closest point.
FESurfaceElement* FEClosestPointProjection::Project(int n, vec3d& q, vec2d& r)
{
	// get the mesh
	FEMesh& mesh = *m_surf.GetMesh();

	// get the node's position
	vec3d x = mesh.Node(n).m_rt;
	
	// see if we need to initialize the NQ structure
//	if (binit_nq) m_SNQ.Init();
	
	// let's find the closest node
//	int mn = m_SNQ.Find(x);

	// do it the hard way
	int mn = -1;
	double d0;
	int N = m_surf.Nodes();
	for (int i=0; i<N; ++i)
	{
		if (m_surf.NodeIndex(i) != n)
		{
			vec3d r = m_surf.Node(i).m_rt;
			double d = (r - x)*(r - x);
			if ((mn == -1) || (d < d0))
			{
				d0 = d;
				mn = i;
			}
		}
	}
	
	// mn is a local index, so get the global node number too
	int m = m_surf.NodeIndex(mn);
	
	// get the nodal position
	vec3d r0 = mesh.Node(m).m_rt;

	// check the distance
	double D = (x - r0).norm();
	if ((m_rad > 0) && (D > m_rad)) return 0;

	// The node cannot be part of the net of the closest point
	int nval = m_NEL.Valence(mn);
	FEElement** pe = m_NEL.ElementList(mn);
	for (int j=0; j<nval; ++j)
	{
		// get the element
		FESurfaceElement& el = static_cast<FESurfaceElement&> (*pe[j]);
		// make sure the node is not contained in this surface element
		if (el.HasNode(n)) return 0;
	}
	
	// now that we found the closest node, lets see if we can find 
	// the best element
	for (int j=0; j<nval; ++j)
	{
		// get the element
		FESurfaceElement& el = static_cast<FESurfaceElement&> (*pe[j]);
		int N = el.Nodes();

		// project the node on the element
		r[0] = 0;
		r[1] = 0;
		q = m_surf.ProjectToSurface(el, x, r[0], r[1]);
		if (m_surf.IsInsideElement(el, r[0], r[1], m_tol)) return &el;
	}
	
	// If we get here, we did not find a facet.
	// There are a couple of reasons why the search has failed:
	// -1. the point cannot be projected onto the surface. For contact this implies the node is not in contact.
	// -2. the projection falls outside the set of elements surrounding the closest point.
	// -3. the projection falls on an edge of two faces whos normals are pointing away.
	// -4. the closest node is in fact the closest point and no closer projection on face or edge can be found
	//
	// TODO: I am not sure yet how to distinguish these cases. One way might be to project the point onto
	//       the edges and see if there is an edge that has a projection that is closer than the closest
	//       node. I am not sure yet if this is a sufficient argument though.
	return 0;
}
