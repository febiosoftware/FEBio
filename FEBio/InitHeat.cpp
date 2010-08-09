#include "stdafx.h"
#include "fem.h"

//-----------------------------------------------------------------------------
//! Initialize heat-conduction data

bool FEM::InitHeat()
{
	int i;

	if (m_pStep->m_nModule != FE_HEAT)
	{
		// if there is no heat-conduction
		// we set all temperature degrees of freedoms as fixed
		// just to make sure they do not get assigned an equation number
		for (i=0; i<m_mesh.Nodes(); ++i) m_mesh.Node(i).m_ID[10] = -1;

		// also remove prescribed temperatures
		for (i=0; i<(int) m_DC.size(); ++i)
		{
			int& node = m_DC[i].node;
			int& bc   = m_DC[i].bc;

			if (bc == 10) bc = -1;
		}

		// let's go back
		return true;
	}

	// for now we fix all non-temperature degrees of freedom
	for (i=0; i<m_mesh.Nodes(); ++i)
	{
		FENode& node = m_mesh.Node(i);
		for (int j=0; j<MAX_NDOFS; ++j)
			if (j != 10) node.m_ID[j] = -1;
	}

	// all done
	return true;
}
