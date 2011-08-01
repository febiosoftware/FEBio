#include "stdafx.h"
#include "fem.h"
#include "FEPoroElastic.h"
#include "FECore/FEMaterial.h"
#include "FEBiphasic.h"
#include "FEElasticShellDomain.h"
#include "FEElasticSolidDomain.h"
#include "log.h"

//-----------------------------------------------------------------------------
//! Initialize solute-poroelastic data

bool FEM::InitPoroSolute()
{
	int i, j, nd;

	// make sure this is the poro-solute module
	bool bporo = (m_pStep->m_nModule == FE_POROELASTIC) || (m_pStep->m_nModule == FE_POROSOLUTE);
	bool bsolu = (m_pStep->m_nModule == FE_POROSOLUTE);
	
	if (!bporo)
	{
		// if there is no poroelasticity
		// we set all pressure degrees of freedoms as fixed
		// just to make sure they do not get assigned an equation number
		for (i=0; i<m_mesh.Nodes(); ++i) m_mesh.Node(i).m_ID[DOF_P] = -1;
		
		// also remove prescribed pressures
		for (i=0; i<(int) m_DC.size(); ++i)
		{
			int& bc   = m_DC[i]->bc;
			if (bc == DOF_P) bc = -1;
		}
	}
	
	if (!bsolu)
	{
		// if there is no solute
		// we set all concentration degrees of freedoms as fixed
		// just to make sure they do not get assigned an equation number
		for (i=0; i<m_mesh.Nodes(); ++i) m_mesh.Node(i).m_ID[DOF_C] = -1;
		
		// also remove prescribed concentrations
		for (i=0; i<(int) m_DC.size(); ++i)
		{
			int& bc   = m_DC[i]->bc;
			if (bc == DOF_C) bc = -1;
		}
	}
	
	if ((!bporo) && (!bsolu))
		// let's go back
		return true;
	
	// see if we are using the symmetric version or not
	if (m_bsym_poro == false) 
	{
		SetSymmetryFlag(false);
		
		// make sure we are using full-Newton
//		if (m_pStep->m_psolver->m_bfgs.m_maxups != 0)
//		{
//			m_pStep->m_psolver->m_bfgs.m_maxups = 0;
//			clog.printbox("WARNING", "The non-symmetric solver algorithm does not work with BFGS yet.\nThe full-Newton method will be used instead.");
//		}
	}
	
	if (bporo)
	{
		// fix all pressure dofs that are not used
		// that is, that are not part of a poroelastic or solute-poroelastic element
		// this is done in three steps
		// step 1. mark all poroelastic nodes
		for (nd = 0; nd<m_mesh.Domains(); ++nd)
		{
			FEElasticSolidDomain* pbd = dynamic_cast<FEElasticSolidDomain*>(&m_mesh.Domain(nd));
			if (pbd)
			{
				for (i=0; i<pbd->Elements(); ++i)
				{
					FESolidElement& el = pbd->Element(i);
					FEPoroElastic* pm = dynamic_cast<FEPoroElastic*>(GetMaterial(el.GetMatID()));
					FEBiphasic* bm = dynamic_cast<FEBiphasic*>(GetMaterial(el.GetMatID()));
					FEBiphasicSolute* bsm = dynamic_cast<FEBiphasicSolute*>(GetMaterial(el.GetMatID()));
					if (pm || bm || bsm)
					{
						int N = el.Nodes();
						int* n = &el.m_node[0];
						for (j=0; j<N; ++j) 
							if (m_mesh.Node(n[j]).m_ID[DOF_P] == 0) m_mesh.Node(n[j]).m_ID[DOF_P] = 1;
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
					FEBiphasic* bm = dynamic_cast<FEBiphasic*>(GetMaterial(el.GetMatID()));
					FEBiphasicSolute* bsm = dynamic_cast<FEBiphasicSolute*>(GetMaterial(el.GetMatID()));
					if (pm || bm || bsm)
					{
						int N = el.Nodes();
						int* n = &el.m_node[0];
						for (j=0; j<N; ++j) 
							if (m_mesh.Node(n[j]).m_ID[DOF_P] == 0) m_mesh.Node(n[j]).m_ID[DOF_P] = 1;
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
					FEBiphasic* bm = dynamic_cast<FEBiphasic*>(GetMaterial(el.GetMatID()));
					FEBiphasicSolute* bsm = dynamic_cast<FEBiphasicSolute*>(GetMaterial(el.GetMatID()));
					if ((pm == 0) && (bm == 0) && (bsm == 0))
					{
						int N = el.Nodes();
						int* n = &el.m_node[0];
						for (j=0; j<N; ++j)
							if (m_mesh.Node(n[j]).m_ID[DOF_P] != 1) m_mesh.Node(n[j]).m_ID[DOF_P] = -1;
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
					FEBiphasic* bm = dynamic_cast<FEBiphasic*>(GetMaterial(el.GetMatID()));
					FEBiphasicSolute* bsm = dynamic_cast<FEBiphasicSolute*>(GetMaterial(el.GetMatID()));
					if ((pm == 0) && (bm == 0) && (bsm == 0))
					{
						int N = el.Nodes();
						int* n = &el.m_node[0];
						for (j=0; j<N; ++j)
							if (m_mesh.Node(n[j]).m_ID[DOF_P] != 1) m_mesh.Node(n[j]).m_ID[DOF_P] = -1;
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
					FEBiphasic* bm = dynamic_cast<FEBiphasic*>(GetMaterial(el.GetMatID()));
					FEBiphasicSolute* bsm = dynamic_cast<FEBiphasicSolute*>(GetMaterial(el.GetMatID()));
					if (pm || bm || bsm)
					{
						int N = el.Nodes();
						int* n = &el.m_node[0];
						for (j=0; j<N; ++j)
							if (m_mesh.Node(n[j]).m_ID[DOF_P] == 1) m_mesh.Node(n[j]).m_ID[DOF_P] = 0;
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
					FEBiphasic* bm = dynamic_cast<FEBiphasic*>(GetMaterial(el.GetMatID()));
					FEBiphasicSolute* bsm = dynamic_cast<FEBiphasicSolute*>(GetMaterial(el.GetMatID()));
					if (pm || bm || bsm)
					{
						int N = el.Nodes();
						int* n = &el.m_node[0];
						for (j=0; j<N; ++j)
							if (m_mesh.Node(n[j]).m_ID[DOF_P] == 1) m_mesh.Node(n[j]).m_ID[DOF_P] = 0;
					}
				}
			}
		}
	}
	
	if (bsolu)
	{
		// fix all concentration dofs that are not used
		// that is, that are not part of a solute-solid element
		// this is done in three steps
		// step 1. mark all solute-solid nodes
		for (nd = 0; nd<m_mesh.Domains(); ++nd)
		{
			FEElasticSolidDomain* pbd = dynamic_cast<FEElasticSolidDomain*>(&m_mesh.Domain(nd));
			if (pbd)
			{
				for (i=0; i<pbd->Elements(); ++i)
				{
					FESolidElement& el = pbd->Element(i);
					FEBiphasicSolute* bsm = dynamic_cast<FEBiphasicSolute*>(GetMaterial(el.GetMatID()));
					if (bsm)
					{
						int N = el.Nodes();
						int* n = &el.m_node[0];
						for (j=0; j<N; ++j) 
							if (m_mesh.Node(n[j]).m_ID[DOF_C] == 0) m_mesh.Node(n[j]).m_ID[DOF_C] = 1;
					}
				}
			}
			
			FEElasticShellDomain* psd = dynamic_cast<FEElasticShellDomain*>(&m_mesh.Domain(nd));
			if (psd)
			{
				for (i=0; i<psd->Elements(); ++i)
				{
					FEShellElement& el = psd->Element(i);
					FEBiphasicSolute* bsm = dynamic_cast<FEBiphasicSolute*>(GetMaterial(el.GetMatID()));
					if (bsm)
					{
						int N = el.Nodes();
						int* n = &el.m_node[0];
						for (j=0; j<N; ++j) 
							if (m_mesh.Node(n[j]).m_ID[DOF_C] == 0) m_mesh.Node(n[j]).m_ID[DOF_C] = 1;
					}
				}
			}
		}
		
		// step 2. fix concentration dofs of all unmarked nodes
		for (nd = 0; nd<m_mesh.Domains(); ++nd)
		{
			FEElasticSolidDomain* pbd = dynamic_cast<FEElasticSolidDomain*>(&m_mesh.Domain(nd));
			if (pbd)
			{
				for (i=0; i<pbd->Elements(); ++i)
				{
					FESolidElement& el = pbd->Element(i);
					FEBiphasicSolute* bsm = dynamic_cast<FEBiphasicSolute*>(GetMaterial(el.GetMatID()));
					if (bsm == 0)
					{
						int N = el.Nodes();
						int* n = &el.m_node[0];
						for (j=0; j<N; ++j)
							if (m_mesh.Node(n[j]).m_ID[DOF_C] != 1) m_mesh.Node(n[j]).m_ID[DOF_C] = -1;
					}
				}
			}
			
			FEElasticShellDomain* psd = dynamic_cast<FEElasticShellDomain*>(&m_mesh.Domain(nd));
			if (psd)
			{
				for (i=0; i<psd->Elements(); ++i)
				{
					FEShellElement& el = psd->Element(i);
					FEBiphasicSolute* bsm = dynamic_cast<FEBiphasicSolute*>(GetMaterial(el.GetMatID()));
					if (bsm == 0)
					{
						int N = el.Nodes();
						int* n = &el.m_node[0];
						for (j=0; j<N; ++j)
							if (m_mesh.Node(n[j]).m_ID[DOF_C] != 1) m_mesh.Node(n[j]).m_ID[DOF_C] = -1;
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
					FEBiphasicSolute* bsm = dynamic_cast<FEBiphasicSolute*>(GetMaterial(el.GetMatID()));
					if (bsm)
					{
						int N = el.Nodes();
						int* n = &el.m_node[0];
						for (j=0; j<N; ++j)
							if (m_mesh.Node(n[j]).m_ID[DOF_C] == 1) m_mesh.Node(n[j]).m_ID[DOF_C] = 0;
					}
				}
			}
			
			
			FEElasticShellDomain* psd = dynamic_cast<FEElasticShellDomain*>(&m_mesh.Domain(nd));
			if (psd)
			{
				for (i=0; i<psd->Elements(); ++i)
				{
					FEShellElement& el = psd->Element(i);
					FEBiphasicSolute* bsm = dynamic_cast<FEBiphasicSolute*>(GetMaterial(el.GetMatID()));
					if (bsm)
					{
						int N = el.Nodes();
						int* n = &el.m_node[0];
						for (j=0; j<N; ++j)
							if (m_mesh.Node(n[j]).m_ID[DOF_C] == 1) m_mesh.Node(n[j]).m_ID[DOF_C] = 0;
					}
				}
			}
		}
	}
	
	return true;
}
