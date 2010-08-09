#include "stdafx.h"
#include "fem.h"
#include "FEPoroElastic.h"
#include "log.h"

//-----------------------------------------------------------------------------
//! Initialize poro-elastic data

bool FEM::InitPoro()
{
	int i, j, nd;

	// see if there are any poro-elastic materials present
	for (i=0; i<Materials(); ++i)
		if (dynamic_cast<FEPoroElastic*>(&m_MAT[i]))
		{
			m_pStep->m_nModule = FE_POROELASTIC;
			break;
		}

	if (m_pStep->m_nModule != FE_POROELASTIC)
	{
		// if there is no poroelasticity
		// we set all pressure degrees of freedoms as fixed
		// just to make sure they do not get assigned an equation number
		for (i=0; i<m_mesh.Nodes(); ++i) m_mesh.Node(i).m_ID[6] = -1;

		// also remove prescribed pressures
		for (i=0; i<(int) m_DC.size(); ++i)
		{
			int& node = m_DC[i].node;
			int& bc   = m_DC[i].bc;

			if (bc == 6) bc = -1;
		}

		// let's go back
		return true;
	}

	// get the logfile
	Logfile& log = GetLogfile();

	// see if we are the symmetric version or not
	if (m_bsym_poro == false) 
	{
		SetSymmetryFlag(false);
	
		// make sure we are using full-Newton
		if (m_pStep->m_psolver->m_bfgs.m_maxups != 0)
		{
			m_pStep->m_psolver->m_bfgs.m_maxups = 0;
			log.printbox("WARNING", "The non-symmetric biphasic algorithm does not work with BFGS yet.\nThe full-Newton method will be used instead.");
		}
	}

	// fix all pressure dofs that are not used
	// that is, that are not part of a poro-elastic element
	// this is done in three steps
	// step 1. mark all poro-elastic nodes
	for (nd = 0; nd<m_mesh.Domains(); ++nd)
	{
		FEElasticSolidDomain* pbd = dynamic_cast<FEElasticSolidDomain*>(&m_mesh.Domain(nd));
		if (pbd)
		{
			for (i=0; i<pbd->Elements(); ++i)
			{
				FESolidElement& el = pbd->Element(i);
				FEPoroElastic* pm = dynamic_cast<FEPoroElastic*>(GetMaterial(el.GetMatID()));
				if (pm)
				{
					int N = el.Nodes();
					int* n = &el.m_node[0];
					for (j=0; j<N; ++j) 
						if (m_mesh.Node(n[j]).m_ID[6] == 0) m_mesh.Node(n[j]).m_ID[6] = 1;
				}
			}
		}

		FEElasticShellDomain* psd = dynamic_cast<FEElasticShellDomain*>(&m_mesh.Domain(nd));
		if (psd)
		{
			for (i=0; i<psd->Elements(); ++i)
			{
				FEShellElement& el = psd->Element(i);
				FEPoroElastic* pm = dynamic_cast<FEPoroElastic*>(GetMaterial(el.GetMatID()));
				if (pm)
				{
					int N = el.Nodes();
					int* n = &el.m_node[0];
					for (j=0; j<N; ++j) 
						if (m_mesh.Node(n[j]).m_ID[6] == 0) m_mesh.Node(n[j]).m_ID[6] = 1;
				}
			}
		}
	}

	// step 2. fix pressure dofs of all unmarked nodes
	for (nd = 0; nd<m_mesh.Domains(); ++nd)
	{
		FEElasticSolidDomain* pbd = dynamic_cast<FEElasticSolidDomain*>(&m_mesh.Domain(nd));
		if (pbd)
		{
			for (i=0; i<pbd->Elements(); ++i)
			{
				FESolidElement& el = pbd->Element(i);
				FEPoroElastic* pm = dynamic_cast<FEPoroElastic*>(GetMaterial(el.GetMatID()));
				if (pm == 0)
				{
					int N = el.Nodes();
					int* n = &el.m_node[0];
					for (j=0; j<N; ++j)
						if (m_mesh.Node(n[j]).m_ID[6] != 1) m_mesh.Node(n[j]).m_ID[6] = -1;
				}
			}
		}

		FEElasticShellDomain* psd = dynamic_cast<FEElasticShellDomain*>(&m_mesh.Domain(nd));
		if (psd)
		{
			for (i=0; i<psd->Elements(); ++i)
			{
				FEShellElement& el = psd->Element(i);
				FEPoroElastic* pm = dynamic_cast<FEPoroElastic*>(GetMaterial(el.GetMatID()));
				if (pm == 0)
				{
					int N = el.Nodes();
					int* n = &el.m_node[0];
					for (j=0; j<N; ++j)
						if (m_mesh.Node(n[j]).m_ID[6] != 1) m_mesh.Node(n[j]).m_ID[6] = -1;
				}
			}
		}
	}

	// step 3. free all marked dofs
	for (nd = 0; nd<m_mesh.Domains(); ++nd)
	{
		FEElasticSolidDomain* pbd = dynamic_cast<FEElasticSolidDomain*>(&m_mesh.Domain(nd));
		if (pbd)
		{
			for (i=0; i<pbd->Elements(); ++i)
			{
				FESolidElement& el = pbd->Element(i);
				FEPoroElastic* pm = dynamic_cast<FEPoroElastic*>(GetMaterial(el.GetMatID()));
				if (pm)
				{
					int N = el.Nodes();
					int* n = &el.m_node[0];
					for (j=0; j<N; ++j)
						if (m_mesh.Node(n[j]).m_ID[6] == 1) m_mesh.Node(n[j]).m_ID[6] = 0;
				}
			}
		}


		FEElasticShellDomain* psd = dynamic_cast<FEElasticShellDomain*>(&m_mesh.Domain(nd));
		if (psd)
		{
			for (i=0; i<psd->Elements(); ++i)
			{
				FEShellElement& el = psd->Element(i);
				FEPoroElastic* pm = dynamic_cast<FEPoroElastic*>(GetMaterial(el.GetMatID()));
				if (pm)
				{
					int N = el.Nodes();
					int* n = &el.m_node[0];
					for (j=0; j<N; ++j)
						if (m_mesh.Node(n[j]).m_ID[6] == 1) m_mesh.Node(n[j]).m_ID[6] = 0;
				}
			}
		}
	}

	// determined the nr of pressure equations
	m_npeq = 0;
	for (i=0; i<m_mesh.Nodes(); ++i)
	{
		FENode& n = m_mesh.Node(i);
		if (n.m_ID[6] != -1) m_npeq++;
	}

	return true;
}
