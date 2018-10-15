// FEBox.cpp: implementation of the FEBox class.
//
//////////////////////////////////////////////////////////////////////

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

void FEBoxMesh::Create(int nx, int ny, int nz, vec3d r0, vec3d r1, int nhex)
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

				// set rigid body id
				node.m_rid = -1;
			}
		}
	}

	// create the elements
	int *en;
	n = 0;
	FEElasticSolidDomain* pbd = new FEElasticSolidDomain(fem);
	pbd->Create(elems, nhex);
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
