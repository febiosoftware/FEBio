#include "stdafx.h"
#include "FESolidSolver.h"
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

	// Get the logfile
	Logfile& logf = GetLogfile();

	// Do rigid joint augmentations
	if (m_fem.m_nrj)
	{
		// loop over all rigid joints
		for (int i=0; i<m_fem.m_nrj; ++i) bconv = m_fem.m_RJ[i].Augment() && bconv;
	}

	// Do contact augmentations
	if (m_fem.m_bcontact)
	{
		// loop over all contact interfaces
		for (int i=0; i<m_fem.ContactInterfaces(); ++i) bconv = m_fem.m_CI[i].Augment(m_naug) && bconv;
	}

	// do linear constraint augmentations
	if (m_fem.m_LCSet.size())
	{
		int n = m_fem.m_LCSet.size();
		list<FELinearConstraintSet*>::iterator im = m_fem.m_LCSet.begin();
		for (int i=0; i<n; ++i, ++im) bconv = (*im)->Augment(m_naug) && bconv;
	}

	// do incompressibility multipliers
	for (int i=0; i<m_fem.Materials(); ++i)
	{
		FEIncompressibleMaterial* pmi = dynamic_cast<FEIncompressibleMaterial*>(m_fem.GetMaterial(i));
		if (pmi)
		{
			if (pmi->m_blaugon)
			{
				int n;
				double normL0 = 0, normL1 = 0, L0, L1;
				double k = pmi->m_K;
				FEMesh& mesh = m_fem.m_mesh;

				// do solid elements
				FESolidDomain& bd = mesh.SolidDomain();
				for (n=0; n<bd.Elements(); ++n)
				{
					FESolidElement& el = bd.Element(n);

					if (el.GetMatID() == i)
					{
						L0 = el.m_Lk;
						normL0 += L0*L0;

						L1 = L0 + k*pmi->h(el.m_eJ);
						normL1 += L1*L1;
					}
				}

				// do shell elements
				FEShellDomain& sd = mesh.ShellDomain();
				for (n=0; n<sd.Elements(); ++n)
				{
					FEShellElement& el = sd.Element(n);

					if (el.GetMatID() == i)
					{
						L0 = el.m_Lk;
						normL0 += L0*L0;

						L1 = L0 + k*pmi->h(el.m_eJ);
						normL1 += L1*L1;
					}
				}
				normL0 = sqrt(normL0);
				normL1 = sqrt(normL1);

				// check convergence
				double pctn = 0;
				if (fabs(normL1) > 1e-10) pctn = fabs((normL1 - normL0)/normL1);

				logf.printf(" material %d\n", i+1);
				logf.printf("                        CURRENT         CHANGE        REQUIRED\n");
				logf.printf("   pressure norm : %15le%15le%15le\n", normL1, pctn, pmi->m_atol);

				if (pctn >= pmi->m_atol)
				{
					bconv = false;
					FESolidDomain& bd = mesh.SolidDomain();
					for (n=0; n<bd.Elements(); ++n)
					{
						FESolidElement& el = bd.Element(n);
						if (el.GetMatID() == i) 
						{
							double hi = pmi->h(el.m_eJ);
							el.m_Lk += k*pmi->h(el.m_eJ);
							el.m_ep = el.m_Lk*pmi->hp(el.m_eJ) + k*log(el.m_eJ)/el.m_eJ;
						}
					}

					FEShellDomain& sd = mesh.ShellDomain();
					for (n=0; n<sd.Elements(); ++n)
					{
						FEShellElement& el = sd.Element(n);
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

	return bconv;
}
