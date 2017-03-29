#include "stdafx.h"
#include "FEMergedConstraint.h"
#include <FECore/FEModel.h>
#include <FECore/FELinearConstraintManager.h>
#include <FECore/FELinearConstraint.h>

FEMergedConstraint::FEMergedConstraint(FEModel& fem) : m_fem(fem)
{
}

FEMergedConstraint::~FEMergedConstraint()
{
}

bool FEMergedConstraint::Merge(FEFacetSet* surf1, FEFacetSet* surf2, const vector<int>& dofList)
{
	// extract the nodes from the surfaces
	FENodeSet set1 = surf1->GetNodeSet();
	FENodeSet set2 = surf2->GetNodeSet();

	// find for each node on surface1 a corresponding node on surface 2 within tolerance
	// First, make sure that set2 is larger than set1
	if (set1.size() > set2.size()) return false;
	int N1 = set1.size();
	int N2 = set2.size();

	// make sure there is something to do
	if (N1 == 0) return true;
	if (dofList.size() == 0) return true;

	// alright, let's get going
	vector<int> tag(N1, -1);
	for (int i=0; i<N1; ++i)
	{
		// get the node position
		vec3d ri = set1.Node(i)->m_rt;

		// find the closest node
		int n = 0;
		double Dmin = (set2.Node(0)->m_rt - ri).norm2();
		for (int j=1; j<N2; ++j)
		{
			vec3d rj = set2.Node(j)->m_rt;
			double D2 = (ri - rj).norm2();
			if (D2 < Dmin)
			{
				Dmin = D2;
				n = j;
			}
		}

		// since this interface type assumes that the nodes match identically,
		// we check that the min distance is indeed very small
		if (Dmin > 1e-9) { 
			return false;
		}

		// store this node
		tag[i] = n;
	}

	// next, create the linear constraints 
	// get the linear constraint manager
	int ndofs = (int) dofList.size();
	FELinearConstraintManager& LCM = m_fem.GetLinearConstraintManager();
	for (int i=0; i<N1; ++i)
	{
		for (int j=0; j<ndofs; ++j)
		{
			int dof = dofList[j];
			FELinearConstraint lc(&m_fem);
			lc.SetMasterDOF(dof, set1[i]);
			lc.AddSlaveDof(dof, set2[tag[i]], 1.0);

			LCM.AddLinearConstraint(lc);
		}
	}

	// all done
	return true;
}
