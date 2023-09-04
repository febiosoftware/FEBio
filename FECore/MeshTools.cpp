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
#include "MeshTools.h"
#include "FEBoundingBox.h"
#include "FEMesh.h"

// Build the bounding box of a surface
FEBoundingBox CalculateBoundingBox(FESurface* ps)
{
	FEBoundingBox box;
	if (ps == nullptr) return box;

	for (int i = 0; i < ps->Nodes(); ++i)
	{
		FENode& node = ps->Node(i);
		if (i == 0) { box = FEBoundingBox(node.m_rt); }
		else box.add(node.m_rt);
	}

	return box;
}

// Calculate the intersected edges of a domain and an immersed boundary
FEEdgeList FindIntersectedEdges(FEDomain* dom, FESurface* ps)
{
	// we'll need the mesh
	FEMesh* pm = (dom ? dom->GetMesh() : nullptr); assert(pm);

	// the intersected edge list that we'll return 
	FEEdgeList IEL(pm);
	if (pm == nullptr) return IEL;

	// Build the complete edge list of the domain
	FEEdgeList EL; 
	if (EL.Create(dom) == false) { assert(false); return IEL; }

	// calculate a bounding box of the surface so we can quickly
	// eliminate edges that definitely won't intersect the surface. 
	FEBoundingBox box = CalculateBoundingBox(ps);

	// inflate it a little to avoid any numerical rounding issues
	double dR = 1e-5*box.radius(); 
	box.inflate(dR, dR, dR);

	// loop over all the edges
	for (int i = 0; i < EL.Edges(); ++i)
	{
		const FEEdgeList::EDGE& e = EL[i];

		// get the two node positions
		vec3d r0 = dom->Node(e.node[0]).m_rt;
		vec3d r1 = dom->Node(e.node[1]).m_rt;

		int n0 = dom->NodeIndex(e.node[0]);
		int n1 = dom->NodeIndex(e.node[1]);

		// do a quick test to see of this edge intersects the bounding box
		// (There might be cases where neither point is in the box, but the edge
		//  still intersects it. May need to expand this test.)
		if (box.IsInside(r0) || box.IsInside(r1))
		{
			vec3d n = r1 - r0;

			// see if this edge intersects with the surface
			int NF = ps->Elements();
			for (int j = 0; j < NF; ++j)
			{
				FESurfaceElement& el = ps->Element(j);

				// find the intersection with the element
				// This function will only return true if the ray
				// intersects from the positive side, so we need to test twice
				double rs[2] = { 0 }, g(0.0), eps(1e-5);
				if (ps->Intersect(el, r0, n, rs, g, eps))
				{
					// we found an intersection with the ray (r0, n), but
					// does the edge actually intersect the surface? 
					if ((g > 0) && (g < 1))
					{
						IEL.Add(n0, n1, g);
					}
				}
				else if (ps->Intersect(el, r1, -n, rs, g, eps))
				{
					// we found an intersection with the ray (r0, n), but
					// does the edge actually intersect the surface? 
					if ((g > 0) && (g < 1))
					{
						IEL.Add(n1, n0, 1.0 - g);
					}
				}
			}
		}
	}

	// all done
	return IEL;
}