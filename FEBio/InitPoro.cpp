#include "stdafx.h"
#include "fem.h"

//-----------------------------------------------------------------------------
//! Initialize poro-elastic data

bool FEM::InitPoro()
{
	int i, j;

	if (m_pStep->m_itype != FE_STATIC_PORO)
	{
		// if there is no poroelasticity
		// we set all pressure degrees of freedoms as fixed
		// just to make sure they do not get assigned an equation number
		for (i=0; i<m_mesh.Nodes(); ++i) m_mesh.Node(i).m_ID[6] = -1;

		// also remove prescribed pressures
		for (i=0; i<m_DC.size(); ++i)
		{
			int& node = m_DC[i].node;
			int& bc   = m_DC[i].bc;

			if (bc == 6) bc = -1;
		}

		// let's go back
		return true;
	}

	// fix all pressure dofs that are not used
	// that is, that are not part of a poro-elastic element
	// this is done in three steps
	// step 1. mark all poro-elastic nodes
	for (i=0; i<m_mesh.SolidElements(); ++i)
	{
		FESolidElement& el = m_mesh.SolidElement(i);
		FEMaterial* pm = GetMaterial(el.GetMatID());
		if (pm->Type() == FE_PORO_ELASTIC)
		{
			int N = el.Nodes();
			int* n = el.m_node;
			for (j=0; j<N; ++j) 
				if (m_mesh.Node(n[j]).m_ID[6] == 0) m_mesh.Node(n[j]).m_ID[6] = 1;
		}
	}
	for (i=0; i<m_mesh.ShellElements(); ++i)
	{
		FEShellElement& el = m_mesh.ShellElement(i);
		FEMaterial* pm = GetMaterial(el.GetMatID());
		if (pm->Type() == FE_PORO_ELASTIC)
		{
			int N = el.Nodes();
			int* n = el.m_node;
			for (j=0; j<N; ++j) 
				if (m_mesh.Node(n[j]).m_ID[6] == 0) m_mesh.Node(n[j]).m_ID[6] = 1;
		}
	}

	// step 2. fix pressure dofs of all unmarked nodes
	for (i=0; i<m_mesh.SolidElements(); ++i)
	{
		FESolidElement& el = m_mesh.SolidElement(i);
		FEMaterial* pm = GetMaterial(el.GetMatID());
		if (pm->Type() != FE_PORO_ELASTIC)
		{
			int N = el.Nodes();
			int* n = el.m_node;
			for (j=0; j<N; ++j)
				if (m_mesh.Node(n[j]).m_ID[6] != 1) m_mesh.Node(n[j]).m_ID[6] = -1;
		}
	}
	for (i=0; i<m_mesh.ShellElements(); ++i)
	{
		FEShellElement& el = m_mesh.ShellElement(i);
		FEMaterial* pm = GetMaterial(el.GetMatID());
		if (pm->Type() != FE_PORO_ELASTIC)
		{
			int N = el.Nodes();
			int* n = el.m_node;
			for (j=0; j<N; ++j)
				if (m_mesh.Node(n[j]).m_ID[6] != 1) m_mesh.Node(n[j]).m_ID[6] = -1;
		}
	}

	// step 3. free all marked dofs
	for (i=0; i<m_mesh.SolidElements(); ++i)
	{
		FESolidElement& el = m_mesh.SolidElement(i);
		FEMaterial* pm = GetMaterial(el.GetMatID());
		if (pm->Type() == FE_PORO_ELASTIC)
		{
			int N = el.Nodes();
			int* n = el.m_node;
			for (j=0; j<N; ++j)
				if (m_mesh.Node(n[j]).m_ID[6] == 1) m_mesh.Node(n[j]).m_ID[6] = 0;
		}
	}
	for (i=0; i<m_mesh.ShellElements(); ++i)
	{
		FEShellElement& el = m_mesh.ShellElement(i);
		FEMaterial* pm = GetMaterial(el.GetMatID());
		if (pm->Type() == FE_PORO_ELASTIC)
		{
			int N = el.Nodes();
			int* n = el.m_node;
			for (j=0; j<N; ++j)
				if (m_mesh.Node(n[j]).m_ID[6] == 1) m_mesh.Node(n[j]).m_ID[6] = 0;
		}
	}

	return true;
}
