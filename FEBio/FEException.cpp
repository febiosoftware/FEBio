// FEException.cpp: implementation of the FEException class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FEException.h"
#include "fem.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

FEException::FEException()
{

}

FEException::~FEException()
{

}

ZeroDiagonal::ZeroDiagonal(int n, FEM& fem)
{
	// let's find what dof this equation belonges to
	FEMesh& m = fem.m_mesh;

	int i, j, id;
	bool bfound = false;
	for (i=0; i<fem.m_nrb; ++i)
	{
		FERigidBody& rb = fem.m_RB[i];
		for (j=0; j<6; ++j)
		{
			id = rb.m_LM[j];
			if (id < -1) id = -id-2;
			if (id == n)
			{
				sprintf(m_szerr, "Zero diagonal on row %d.\nThis dof belongs to rigid body %d (dof %d)", n, i+1, j+1);
				bfound = true;
				break;
			}
		}
		if (bfound) break;
	}

	if (!bfound)
	{
		for (i=0; i<m.Nodes(); ++i)
		{
			FENode& node = m.Node(i);
			for (j=0; j<MAX_NDOFS; ++j) 
				if (node.m_ID[j] == n)
				{
					sprintf(m_szerr, "Zero diagonal on row %d.\nThis dof belongs to node %d (dof %d)", n, i+1, j+1);
					bfound = true;
					break;
				}
			if (bfound) break;
		}
	}

	if (!bfound)
	{
		sprintf(m_szerr, "Zero diagonal on row %d.", n);
		bfound = true;
	}
}
