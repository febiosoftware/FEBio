#include "stdafx.h"
#include "FESolidSolver.h"
#include "FEUncoupledMaterial.h"
#include "log.h"

///////////////////////////////////////////////////////////////////////////////
// FUNCTION: FESolidSolver::Augment
//  This functions performs the Lagrange augmentations
//  It returns true if all the augmentation have converged, 
//	otherwise it returns false
//
// TODO: There is an inherent problem with this approach. Since
//  Lagrangian multipliers are inherited from previous timesteps
//  they might not be zero in case a node-surface contact breaks. 
//  The node's gap value needs to become negative to a certain value
//  before the Lagr. multipliers dissapears. 
//

bool FESolidSolver::Augment()
{
	// Assume we will pass (can't hurt to be optimistic)
	bool bconv = true;

	// Do rigid joint augmentations
	if (!m_fem.m_RJ.empty())
	{
		// loop over all rigid joints
		for (int i=0; i<(int) m_fem.m_RJ.size(); ++i) bconv = m_fem.m_RJ[i]->Augment() && bconv;
	}

	// Do contact augmentations
	if (m_fem.ContactInterfaces() > 0)
	{
		// loop over all contact interfaces
		for (int i=0; i<m_fem.ContactInterfaces(); ++i) bconv = m_fem.m_CI[i]->Augment(m_naug) && bconv;
	}

	// do linear constraint augmentations
	if (m_fem.m_LCSet.size())
	{
		int n = m_fem.m_LCSet.size();
		list<FELinearConstraintSet*>::iterator im = m_fem.m_LCSet.begin();
		for (int i=0; i<n; ++i, ++im) bconv = (*im)->Augment(m_naug) && bconv;
	}

	// do incompressibility multipliers
	int nd;
	for (int i=0; i<m_fem.Materials(); ++i)
	{
		FEUncoupledMaterial* pmi = dynamic_cast<FEUncoupledMaterial*>(m_fem.GetMaterial(i));
		if (pmi)
		{
			if (pmi->m_blaugon)
			{
				int n;
				double normL0 = 0, normL1 = 0, L0, L1;
				double k = pmi->m_K;
				FEMesh& mesh = m_fem.m_mesh;

				for (nd = 0; nd < mesh.Domains(); ++nd)
				{
					// do solid elements
					FEElasticSolidDomain* pbd = dynamic_cast<FEElasticSolidDomain*>(&mesh.Domain(nd));
					if (pbd)
					{
						for (n=0; n<pbd->Elements(); ++n)
						{
							FESolidElement& el = pbd->Element(n);

							if (el.GetMatID() == i)
							{
								L0 = el.m_Lk;
								normL0 += L0*L0;

								L1 = L0 + k*pmi->h(el.m_eJ);
								normL1 += L1*L1;
							}
						}
					}

					// do shell elements
					FEElasticShellDomain* psd = dynamic_cast<FEElasticShellDomain*>(&mesh.Domain(nd));
					if (psd)
					{
						for (n=0; n<psd->Elements(); ++n)
						{
							FEShellElement& el = psd->Element(n);

							if (el.GetMatID() == i)
							{
								L0 = el.m_Lk;
								normL0 += L0*L0;

								L1 = L0 + k*pmi->h(el.m_eJ);
								normL1 += L1*L1;
							}
						}
					}
				}

				normL0 = sqrt(normL0);
				normL1 = sqrt(normL1);

				// check convergence
				double pctn = 0;
				if (fabs(normL1) > 1e-10) pctn = fabs((normL1 - normL0)/normL1);

				clog.printf(" material %d\n", i+1);
				clog.printf("                        CURRENT         CHANGE        REQUIRED\n");
				clog.printf("   pressure norm : %15le%15le%15le\n", normL1, pctn, pmi->m_atol);

				if (pctn >= pmi->m_atol)
				{
					bconv = false;
					for (nd = 0; nd < mesh.Domains(); ++nd)
					{
						FEElasticSolidDomain* pbd = dynamic_cast<FEElasticSolidDomain*>(&mesh.Domain(nd));
						if (pbd)
						{
							for (n=0; n<pbd->Elements(); ++n)
							{
								FESolidElement& el = pbd->Element(n);
								if (el.GetMatID() == i) 
								{
									double hi = pmi->h(el.m_eJ);
									el.m_Lk += k*pmi->h(el.m_eJ);
									el.m_ep = el.m_Lk*pmi->hp(el.m_eJ) + k*log(el.m_eJ)/el.m_eJ;
								}
							}
						}
	
						FEElasticShellDomain* psd = dynamic_cast<FEElasticShellDomain*>(&mesh.Domain(nd));
						if (psd)
						{
							for (n=0; n<psd->Elements(); ++n)
							{
								FEShellElement& el = psd->Element(n);
								if (el.GetMatID() == i) 
								{
									el.m_Lk += k*pmi->h(el.m_eJ);
									el.m_ep = el.m_Lk*pmi->hp(el.m_eJ) + k*log(el.m_eJ)/el.m_eJ;
								}
							}
						}
					}
				}
			}
		}
	}

	return bconv;
}
