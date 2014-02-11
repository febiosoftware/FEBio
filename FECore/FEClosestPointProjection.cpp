#include "stdafx.h"
#include "FEClosestPointProjection.h"

//-----------------------------------------------------------------------------
// constructor
FEClosestPointProjection::FEClosestPointProjection(FESurface& s) : m_surf(s)
{
	// set default options
	m_tol = 0.01;
	m_rad = 0.0;	// 0 means don't use search radius

	// calculate node-element list
	m_NET.Create(&m_surf);
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
	int m = m_surf.m_node[mn];

	// get the nodal position
	vec3d r0 = mesh.Node(m).m_rt;

	// now that we found the closest master node, lets see if we can find 
	// the best master element
	int nval = m_NET.Valence(mn);
	FEElement** pe = m_NET.ElementList(mn);
	for (int j=0; j<nval; ++j)
	{
		// get the master element
		FESurfaceElement& el = dynamic_cast<FESurfaceElement&> (*pe[j]);
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
		if (m_surf.m_node[i] != n)
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
	int m = m_surf.m_node[mn];
	
	// get the nodal position
	vec3d r0 = mesh.Node(m).m_rt;

	// check the distance
	double D = (x - r0).norm();
	if ((m_rad > 0) && (D > m_rad)) return 0;

	// The node cannot be part of the net of the closest point
	int nval = m_NET.Valence(mn);
	FEElement** pe = m_NET.ElementList(mn);
	for (int j=0; j<nval; ++j)
	{
		// get the master element
		FESurfaceElement& el = dynamic_cast<FESurfaceElement&> (*pe[j]);
		// make sure the node is not contained in this surface element
		if (el.HasNode(n)) return 0;
	}
	
	// now that we found the closest master node, lets see if we can find 
	// the best master element
	for (int j=0; j<nval; ++j)
	{
		// get the master element
		FESurfaceElement& el = dynamic_cast<FESurfaceElement&> (*pe[j]);
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
