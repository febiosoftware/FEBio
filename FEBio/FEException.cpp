// FEException.cpp: implementation of the FEException class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FEException.h"
#include "fem.h"
#include "log.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

FEException::FEException()
{

}

FEException::~FEException()
{

}

ZeroDiagonal::ZeroDiagonal(vector<int>& l, FEM& fem)
{
	// let's find what dof this equation belonges to
	FEMesh& m = fem.m_mesh;

	Logfile& log = GetLogfile();

	int nz = (int) l.size();
	int i, j, id;
	for (i=0; i<fem.m_nrb; ++i)
	{
		FERigidBody& rb = fem.m_RB[i];
		for (j=0; j<6; ++j)
		{
			id = rb.m_LM[j];
			if (id < -1) id = -id-2;

			for (int k=0; k<nz; ++k)
			{
				int n = l[k];
				if (id == n)
				{
					log.printf("Zero diagonal on row %d.\nThis dof belongs to rigid body %d (dof %d)\n", n, i+1, j+1);
				}
			}
		}
	}

	vector<EQUATION> EQT(fem.m_neq);
	for (i=0; i<fem.m_neq; ++i)
	{
		EQUATION& q = EQT[i];
		q.node = -1;
		q.dof = -1;
	}

	for (i=0; i<m.Nodes(); ++i)
	{
		FENode& node = m.Node(i);
		for (j=0; j<MAX_NDOFS; ++j)
		{
			id = node.m_ID[j];
			if (id < -1) id = -id-2;
			if (id != -1)
			{
				assert(id < fem.m_neq);
				EQT[id].node = i;
				EQT[id].dof = j;
			}
		}
	}

	const int NMAX = 128;
	int N = (nz > NMAX? nz : NMAX), n;
	for (i=0; i<NMAX; ++i)
	{
		n = l[i];
		EQUATION& q = EQT[n];
		log.printf("Zero diagonal on row %d. (node %d, dof %d)\n", n+1, q.node+1, q.dof+1);
	}
	if (nz > NMAX) log.printf("(%d out of %d printed)\n", NMAX, nz);

	// print error message
	sprintf(m_szerr, "FATAL ERROR: %d zero(s) found on diagonal.", l.size());
}
