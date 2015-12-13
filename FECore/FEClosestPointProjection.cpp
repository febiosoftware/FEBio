#include "stdafx.h"
#include "FEClosestPointProjection.h"
#include "FEElemElemList.h"

//-----------------------------------------------------------------------------
// constructor
FEClosestPointProjection::FEClosestPointProjection(FESurface& s) : m_surf(s)
{
	// set default options
	m_tol = 0.01;
	m_rad = 0.0;	// 0 means don't use search radius
	m_bspecial = false;

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

	// let's find the closest master node
	int mn = m_SNQ.Find(x);

	// mn is a local index, so get the global node number too
	int m = m_surf.NodeIndex(mn);

	// get the nodal position
	vec3d rm = mesh.Node(m).m_rt;

	// now that we found the closest master node, lets see if we can find 
	// the best master element
	int nval = m_NEL.Valence(mn);
	FEElement** pe = m_NEL.ElementList(mn);
	int* eil = m_NEL.ElementIndexList(mn);
	for (int j=0; j<nval; ++j)
	{
		// get the master element
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
			// get the master element
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
			// get a master element
			FESurfaceElement& el = static_cast<FESurfaceElement&> (*pe[j]);

			// only consider triangles (for now)
			int N = el.Nodes();
			if (N == 3)
			{
				int kmin = -1;
				// we got to find a point closer than the closest point
				double Dmin = (x - rm).norm();
				vec3d qmin;
				for (int k=0; k<3; ++k)
				{
					// we don't allow projection on boundary edges
					// so make sure the element has a neighbor
					if (m_EEL.Neighbor(eil[j], k))
					{
						// get the two edge node indices
						int nk1 = el.m_node[k];
						int nk2 = el.m_node[(k+1)%3];

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
		
		// first we want to make sure that the node does not lie on the boundary
		for (int j=0; j<nval; ++j)
		{
			// get a master element
			FESurfaceElement& el = static_cast<FESurfaceElement&> (*pe[j]);

			// only consider triangles (for now)
			int N = el.Nodes();
			if (N == 3)
			{
				// make sure that the node does not lie on a boundary edge
				for (int k=0; k<3; ++k)
				{
					if (m_EEL.Neighbor(eil[j], k) == 0)
					{
						int n0 = el.m_node[k];
						int n1 = el.m_node[(k+1)%3];
						if ((n0==m)||(n1==m)) return 0;
					}
				}
			}
		}

		// Any element can now be used but just to be safe, we loop over all
		// and make sure the element actually contains the closest point.
		for (int j=0; j<nval; ++j)
		{
			// get a master element
			FESurfaceElement& el = static_cast<FESurfaceElement&> (*pe[j]);

			// only consider triangles (for now)
			int N = el.Nodes();
			if (N == 3)
			{
				// make sure that the node does not lie on a boundary edge
				for (int k=0; k<3; ++k)
				{
					if (m_EEL.Neighbor(eil[j], k) == 0)
					{
						int n0 = el.m_node[k];
						int n1 = el.m_node[(k+1)%3];
						if ((n0==m)||(n1==m)) return 0;
					}
				}

				int* en = &(el.m_node[0]);

				// make sure one of them is our closest point
				if ((en[0] == m) || (en[1] == m) || (en[2] == m))
				{
					// Now find the iso-parametric coordinates
					q = m_surf.ProjectToSurface(el, rm, r[0], r[1]);
					if (m_surf.IsInsideElement(el, r[0], r[1], m_tol)) return &el;
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
	
	// let's find the closest master node
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
		// get the master element
		FESurfaceElement& el = static_cast<FESurfaceElement&> (*pe[j]);
		// make sure the node is not contained in this surface element
		if (el.HasNode(n)) return 0;
	}
	
	// now that we found the closest master node, lets see if we can find 
	// the best master element
	for (int j=0; j<nval; ++j)
	{
		// get the master element
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
