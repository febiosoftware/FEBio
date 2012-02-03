#include "stdafx.h"
#include "fem.h"
#include "FECore/FEMaterial.h"
#include "FEBioLib/FEBiphasic.h"
#include "FEBioLib/FEBiphasicSolute.h"
#include "FEBioLib/FEElasticShellDomain.h"
#include "FEBioLib/FEElasticSolidDomain.h"
#include "FEBioLib/FETriphasic.h"
#include "FEBioLib/log.h"

//-----------------------------------------------------------------------------
//! Initialize solute-poroelastic data.
//! Find all nodes that are not part of a poro-solute domain and fix the 
//! pressure and concentration DOFS. 
//! \todo This function should probably move to the FEAnalysisStep class.
bool FEM::InitPoroSolute()
{
	int i, j, nd;

	// make sure this is the poro-solute module
	bool bporo = (m_pStep->m_nModule == FE_BIPHASIC) || (m_pStep->m_nModule == FE_POROSOLUTE) || (m_pStep->m_nModule == FE_TRIPHASIC);
	bool bsolu = (m_pStep->m_nModule == FE_POROSOLUTE)  || (m_pStep->m_nModule == FE_TRIPHASIC);
	bool btri  = (m_pStep->m_nModule == FE_TRIPHASIC);
	
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
	
	if (!btri)
	{
		// if there is no solute
		// we set all concentration degrees of freedoms as fixed
		// just to make sure they do not get assigned an equation number
		for (i=0; i<m_mesh.Nodes(); ++i) m_mesh.Node(i).m_ID[DOF_C+1] = -1;
		
		// also remove prescribed concentrations
		for (i=0; i<(int) m_DC.size(); ++i)
		{
			int& bc   = m_DC[i]->bc;
			if (bc == int(DOF_C+1)) bc = -1;
		}
	}
	
	if (!bporo)
		// let's go back
		return true;
	
	// see if we are using the symmetric version or not
	// TODO: this really needs to move to the Analysis step.
	if (m_bsym_poro == false) 
	{
		m_pStep->m_psolver->m_bsymm = false;
		
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
		// that is, that are not part of a biphasic, biphasic-solute, or triphasic element
		// this is done in three steps
		// step 1. mark all biphasic nodes
		for (nd = 0; nd<m_mesh.Domains(); ++nd)
		{
			FEElasticSolidDomain* pbd = dynamic_cast<FEElasticSolidDomain*>(&m_mesh.Domain(nd));
			if (pbd)
			{
				FEBiphasic*       bm  = dynamic_cast<FEBiphasic*      >(pbd->GetMaterial());
				FEBiphasicSolute* bsm = dynamic_cast<FEBiphasicSolute*>(pbd->GetMaterial());
				FETriphasic*	  btm = dynamic_cast<FETriphasic*     >(pbd->GetMaterial());
				for (i=0; i<pbd->Elements(); ++i)
				{
					FESolidElement& el = pbd->Element(i);
					if (bm || bsm || btm)
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
				FEBiphasic*       bm  = dynamic_cast<FEBiphasic*      >(psd->GetMaterial());
				FEBiphasicSolute* bsm = dynamic_cast<FEBiphasicSolute*>(psd->GetMaterial());
				FETriphasic*	  btm = dynamic_cast<FETriphasic*     >(psd->GetMaterial());
				for (i=0; i<psd->Elements(); ++i)
				{
					FEShellElement& el = psd->Element(i);
					if (bm || bsm || btm)
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
				FEBiphasic*       bm  = dynamic_cast<FEBiphasic*      >(pbd->GetMaterial());
				FEBiphasicSolute* bsm = dynamic_cast<FEBiphasicSolute*>(pbd->GetMaterial());
				FETriphasic*	  btm = dynamic_cast<FETriphasic*     >(pbd->GetMaterial());
				for (i=0; i<pbd->Elements(); ++i)
				{
					FESolidElement& el = pbd->Element(i);
					if ((bm == 0) && (bsm == 0) && (btm == 0))
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
				FEBiphasic*       bm  = dynamic_cast<FEBiphasic*      >(psd->GetMaterial());
				FEBiphasicSolute* bsm = dynamic_cast<FEBiphasicSolute*>(psd->GetMaterial());
				FETriphasic*	  btm = dynamic_cast<FETriphasic*     >(psd->GetMaterial());
				for (i=0; i<psd->Elements(); ++i)
				{
					FEShellElement& el = psd->Element(i);
					if ((bm == 0) && (bsm == 0) && (btm == 0))
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
				FEBiphasic*       bm  = dynamic_cast<FEBiphasic*      >(pbd->GetMaterial());
				FEBiphasicSolute* bsm = dynamic_cast<FEBiphasicSolute*>(pbd->GetMaterial());
				FETriphasic*	  btm = dynamic_cast<FETriphasic*     >(pbd->GetMaterial());
				for (i=0; i<pbd->Elements(); ++i)
				{
					FESolidElement& el = pbd->Element(i);
					if (bm || bsm || btm)
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
				FEBiphasic*       bm  = dynamic_cast<FEBiphasic*      >(psd->GetMaterial());
				FEBiphasicSolute* bsm = dynamic_cast<FEBiphasicSolute*>(psd->GetMaterial());
				FETriphasic*	  btm = dynamic_cast<FETriphasic*     >(psd->GetMaterial());
				for (i=0; i<psd->Elements(); ++i)
				{
					FEShellElement& el = psd->Element(i);
					if (bm || bsm || btm)
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
		// fix all neutral/cation concentration dofs that are not used
		// that is, that are not part of a solute-solid element
		// this is done in three steps
		// step 1. mark all solute-solid nodes
		for (nd = 0; nd<m_mesh.Domains(); ++nd)
		{
			FEElasticSolidDomain* pbd = dynamic_cast<FEElasticSolidDomain*>(&m_mesh.Domain(nd));
			if (pbd)
			{
				FEBiphasicSolute* bsm = dynamic_cast<FEBiphasicSolute*>(pbd->GetMaterial());
				FETriphasic*	  btm = dynamic_cast<FETriphasic*     >(pbd->GetMaterial());
				for (i=0; i<pbd->Elements(); ++i)
				{
					FESolidElement& el = pbd->Element(i);
					if (bsm || btm)
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
				FEBiphasicSolute* bsm = dynamic_cast<FEBiphasicSolute*>(psd->GetMaterial());
				FETriphasic*	  btm = dynamic_cast<FETriphasic*     >(psd->GetMaterial());
				for (i=0; i<psd->Elements(); ++i)
				{
					FEShellElement& el = psd->Element(i);
					if (bsm || btm)
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
				FEBiphasicSolute* bsm = dynamic_cast<FEBiphasicSolute*>(pbd->GetMaterial());
				FETriphasic*	  btm = dynamic_cast<FETriphasic*     >(pbd->GetMaterial());
				for (i=0; i<pbd->Elements(); ++i)
				{
					FESolidElement& el = pbd->Element(i);
					if ((bsm == 0) && (btm == 0))
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
				FEBiphasicSolute* bsm = dynamic_cast<FEBiphasicSolute*>(psd->GetMaterial());
				FETriphasic*	  btm = dynamic_cast<FETriphasic*     >(psd->GetMaterial());
				for (i=0; i<psd->Elements(); ++i)
				{
					FEShellElement& el = psd->Element(i);
					if ((bsm == 0) && (btm == 0))
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
				FEBiphasicSolute* bsm = dynamic_cast<FEBiphasicSolute*>(pbd->GetMaterial());
				FETriphasic*	  btm = dynamic_cast<FETriphasic*     >(pbd->GetMaterial());
				for (i=0; i<pbd->Elements(); ++i)
				{
					FESolidElement& el = pbd->Element(i);
					if (bsm || btm)
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
				FEBiphasicSolute* bsm = dynamic_cast<FEBiphasicSolute*>(psd->GetMaterial());
				FETriphasic*	  btm = dynamic_cast<FETriphasic*     >(psd->GetMaterial());
				for (i=0; i<psd->Elements(); ++i)
				{
					FEShellElement& el = psd->Element(i);
					if (bsm || btm)
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

	if (btri)
	{
		// fix all anion concentration dofs that are not used
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
					FETriphasic* btm = dynamic_cast<FETriphasic*>(GetMaterial(el.GetMatID()));
					if (btm)
					{
						int N = el.Nodes();
						int* n = &el.m_node[0];
						for (j=0; j<N; ++j) 
							if (m_mesh.Node(n[j]).m_ID[DOF_C+1] == 0) m_mesh.Node(n[j]).m_ID[DOF_C+1] = 1;
					}
				}
			}
			
			FEElasticShellDomain* psd = dynamic_cast<FEElasticShellDomain*>(&m_mesh.Domain(nd));
			if (psd)
			{
				for (i=0; i<psd->Elements(); ++i)
				{
					FEShellElement& el = psd->Element(i);
					FETriphasic* btm = dynamic_cast<FETriphasic*>(GetMaterial(el.GetMatID()));
					if (btm)
					{
						int N = el.Nodes();
						int* n = &el.m_node[0];
						for (j=0; j<N; ++j) 
							if (m_mesh.Node(n[j]).m_ID[DOF_C+1] == 0) m_mesh.Node(n[j]).m_ID[DOF_C+1] = 1;
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
					FETriphasic* btm = dynamic_cast<FETriphasic*>(GetMaterial(el.GetMatID()));
					if (btm == 0)
					{
						int N = el.Nodes();
						int* n = &el.m_node[0];
						for (j=0; j<N; ++j)
							if (m_mesh.Node(n[j]).m_ID[DOF_C+1] != 1) m_mesh.Node(n[j]).m_ID[DOF_C+1] = -1;
					}
				}
			}
			
			FEElasticShellDomain* psd = dynamic_cast<FEElasticShellDomain*>(&m_mesh.Domain(nd));
			if (psd)
			{
				for (i=0; i<psd->Elements(); ++i)
				{
					FEShellElement& el = psd->Element(i);
					FETriphasic* btm = dynamic_cast<FETriphasic*>(GetMaterial(el.GetMatID()));
					if (btm == 0)
					{
						int N = el.Nodes();
						int* n = &el.m_node[0];
						for (j=0; j<N; ++j)
							if (m_mesh.Node(n[j]).m_ID[DOF_C+1] != 1) m_mesh.Node(n[j]).m_ID[DOF_C+1] = -1;
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
					FETriphasic* btm = dynamic_cast<FETriphasic*>(GetMaterial(el.GetMatID()));
					if (btm)
					{
						int N = el.Nodes();
						int* n = &el.m_node[0];
						for (j=0; j<N; ++j)
							if (m_mesh.Node(n[j]).m_ID[DOF_C+1] == 1) m_mesh.Node(n[j]).m_ID[DOF_C+1] = 0;
					}
				}
			}
			
			
			FEElasticShellDomain* psd = dynamic_cast<FEElasticShellDomain*>(&m_mesh.Domain(nd));
			if (psd)
			{
				for (i=0; i<psd->Elements(); ++i)
				{
					FEShellElement& el = psd->Element(i);
					FETriphasic* btm = dynamic_cast<FETriphasic*>(GetMaterial(el.GetMatID()));
					if (btm)
					{
						int N = el.Nodes();
						int* n = &el.m_node[0];
						for (j=0; j<N; ++j)
							if (m_mesh.Node(n[j]).m_ID[DOF_C+1] == 1) m_mesh.Node(n[j]).m_ID[DOF_C+1] = 0;
					}
				}
			}
		}
	}
	
	return true;
}
