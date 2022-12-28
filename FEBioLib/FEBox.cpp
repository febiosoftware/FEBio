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
#include "FEBox.h"
#include "FEBioMech/FEElasticSolidDomain.h"
#include <FECore/FEModel.h>

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

FEBoxMesh::FEBoxMesh(FEModel* fem) : FEMesh(fem)
{

}

FEBoxMesh::~FEBoxMesh()
{

}

void FEBoxMesh::Create(int nx, int ny, int nz, vec3d r0, vec3d r1, FE_Element_Type nhex)
{
	int i, j, k, n;

	// make sure the parameters make sense
	assert((nx > 0) && (ny > 0) && (nz > 0));

	// count items
	int nodes = (nx+1)*(ny+1)*(nz+1);
	int elems = nx*ny*nz;

	// allocate data
	FEMesh::CreateNodes(nodes);

	FEModel* fem = GetFEModel();
	int MAX_DOFS = fem->GetDOFS().GetTotalDOFS();
	FEMesh::SetDOFS(MAX_DOFS);

	// create the nodes
	double x, y, z;
	n = 0;
	for (i=0; i<=nx; ++i)
	{
		x = r0.x + ((r1.x - r0.x)*i)/nx;
		for (j=0; j<=ny; ++j)
		{
			y = r0.y + ((r1.y - r0.y)*j)/ny;
			for (k=0; k<=nz; ++k, ++n)
			{
				z = r0.z + ((r1.z - r0.z)*k)/nz;

				FENode& node = Node(n);

				node.m_r0 = vec3d(x, y, z);

				node.m_rt = node.m_r0;
			}
		}
	}

	// create the elements
	int *en;
	n = 0;
	FEElasticSolidDomain* pbd = new FEElasticSolidDomain(fem);
	pbd->Create(elems, FEElementLibrary::GetElementSpecFromType(nhex));
	pbd->SetMatID(-1);
	AddDomain(pbd);
	for (i=0; i<nx; ++i)
	{
		for (j=0; j<ny; ++j)
		{
			for (k=0; k<nz; ++k, ++n)
			{
				FESolidElement& el = pbd->Element(n);

				el.SetID(n+1);

				en = &el.m_node[0];

				en[0] = (i  )*(ny+1)*(nz+1) + (j  )*(nz+1) + (k  );
				en[1] = (i+1)*(ny+1)*(nz+1) + (j  )*(nz+1) + (k  );
				en[2] = (i+1)*(ny+1)*(nz+1) + (j+1)*(nz+1) + (k  );
				en[3] = (i  )*(ny+1)*(nz+1) + (j+1)*(nz+1) + (k  );
				en[4] = (i  )*(ny+1)*(nz+1) + (j  )*(nz+1) + (k+1);
				en[5] = (i+1)*(ny+1)*(nz+1) + (j  )*(nz+1) + (k+1);
				en[6] = (i+1)*(ny+1)*(nz+1) + (j+1)*(nz+1) + (k+1);
				en[7] = (i  )*(ny+1)*(nz+1) + (j+1)*(nz+1) + (k+1);
			}
		}
	}
}
