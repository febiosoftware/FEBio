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
#include "FECube.h"
#include "FESurface.h"
#include "FEModel.h"

FECube::FECube() : m_mesh(0)
{
}

FECube::~FECube()
{
	for (int i = 0; i<6; ++i) 
	{
		delete m_surf[i];
		m_surf[i] = 0;
	}
}

// Get a surface
FESurface* FECube::GetSurface(int i)
{
	return m_surf[i];
}

// get the node set of the corner nodes
const FENodeSet& FECube::GetCornerNodes() const
{
	return *m_corners;
}

// get the node set of boundary nodes
const FENodeSet& FECube::GetBoundaryNodes() const
{
	return *m_boundary;
}

// get the mesh of this cube
FEMesh* FECube::GetMesh()
{
	return m_mesh;
}

bool FECube::Build(FEModel* fem)
{
	// make sure we have a mesh
	m_mesh = &fem->GetMesh();

	FEMesh& m = *m_mesh;
	int NN = m.Nodes();
	
	// first, get the outside surface
	FESurface* boundary = m.ElementBoundarySurface();

	// get the boundary node set
	m_boundary = new FENodeSet(fem);
	m_boundary->Add(boundary->GetNodeList());

	// Next, split it up in 6 surfaces
	// We divide the surface by comparing normals to the 6 surface normals of a cube
	vec3d fn[6] = { { 1, 0, 0 }, { -1, 0, 0 }, { 0, 1, 0 }, { 0, -1, 0 }, { 0, 0, 1 }, { 0, 0, -1 } };
	for (int n = 0; n<6; ++n)
	{
		// create the surface
		m_surf[n] = fecore_alloc(FESurface, m_mesh->GetFEModel());

		// get the normal for this face
		vec3d N = fn[n];

		int faces = 0;
		for (int i = 0; i<boundary->Elements(); ++i)
		{
			FESurfaceElement& face = boundary->Element(i);
			vec3d Ni = boundary->SurfaceNormal(face, 0, 0);
			if (Ni*N > 0.9999) faces++;
		}
		m_surf[n]->Create(faces);

		faces = 0;
		for (int i = 0; i<boundary->Elements(); ++i)
		{
			FESurfaceElement& face = boundary->Element(i);
			vec3d Ni = boundary->SurfaceNormal(face, 0, 0);
			if (Ni*N > 0.9999)
			{
				FESurfaceElement& newFace = m_surf[n]->Element(faces++);
				newFace = face;
			}
		}

		m_surf[n]->Init();
	}

	// we also need to find the 8 corner nodes
	vector<int> tag(NN, 0);
	for (int n = 0; n<6; ++n)
	{
		FESurface& sn = *m_surf[n];
		FENodeList ns = sn.GetNodeList();
		for (int i = 0; i<ns.Size(); ++i) tag[ns[i]]++;
	}

	m_corners = new FENodeSet(fem);
	for (int i = 0; i<NN; ++i) if (tag[i] == 3) m_corners->Add(i);
	assert(m_corners->Size() == 8);

	// don't forget to cleanup
	delete boundary;

	return true;
}
