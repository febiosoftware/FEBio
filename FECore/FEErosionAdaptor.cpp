#include "stdafx.h"
#include "FEErosionAdaptor.h"
#include <FECore/FEModel.h>
#include <FECore/FEMesh.h>
#include <FECore/FEDomain.h>
#include <FECore/log.h>
#include <FECore/FELinearConstraintManager.h>
#include <algorithm>

BEGIN_FECORE_CLASS(FEErosionAdaptor, FEMeshAdaptor)
	ADD_PARAMETER(m_maxIters, "max_iters");

	ADD_PROPERTY(m_criterion, "criterion");
END_FECORE_CLASS();

FEErosionAdaptor::FEErosionAdaptor(FEModel* fem) : FEMeshAdaptor(fem)
{
	m_maxIters = -1;
	m_criterion = nullptr;
}

bool FEErosionAdaptor::Apply(int iteration)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	if ((m_maxIters >= 0) && (iteration >= m_maxIters))
	{
		feLog("\tMax iterations reached.");
		return true;
	}

	if (m_criterion == nullptr) return true;

	std::vector<int> selection = m_criterion->GetElementList();
	if (selection.empty())
	{
		feLog("\tNothing to do.\n");
		return true;
	}

	vector<int> elemList(mesh.Elements(), 0);
	for (int i = 0; i < selection.size(); ++i) elemList[selection[i]] = 1;

	int deactiveElems = 0;
	int elem = 0;
	for (int i = 0; i < mesh.Domains(); ++i)
	{
		FEDomain& dom = mesh.Domain(i);
		int NE = dom.Elements();
		for (int j = 0; j < NE; ++j, elem++)
		{
			if (elemList[elem] == 1)
			{
				FEElement& el = dom.ElementRef(j);
				assert(el.isActive());
				el.setInactive();
				deactiveElems++;
			}
		}
	}

	// if any nodes were orphaned, we need to deactivate them as well
	int NN = mesh.Nodes();
	vector<int> tag(NN, 0);
	for (int i = 0; i < mesh.Domains(); ++i)
	{
		FEDomain& dom = mesh.Domain(i);
		int NE = dom.Elements();
		for (int j = 0; j < NE; ++j)
		{
			FEElement& el = dom.ElementRef(j);
			if (el.isActive())
			{
				int neln = el.Nodes();
				for (int n = 0; n < neln; ++n) tag[el.m_node[n]] = 1;
			}
		}
	}

	for (int i = 0; i < NN; ++i)
	{
		FENode& node = mesh.Node(i);
		if (tag[i] == 0)
		{
			node.SetFlags(FENode::EXCLUDE);
			int ndofs = node.dofs();
			for (int j = 0; j < ndofs; ++j)
				node.set_inactive(j);
		}
	}

	// remove any linear constraints of exclude nodes
	FELinearConstraintManager& LCM = fem.GetLinearConstraintManager();
	for (int j = 0; j < LCM.LinearConstraints();)
	{
		FELinearConstraint& lc = LCM.LinearConstraint(j);
		if (mesh.Node(lc.master.node).HasFlags(FENode::EXCLUDE))
		{
			LCM.RemoveLinearConstraint(j);
		}
		else ++j;
	}

	// also remove any linear constraints that have slave excluded nodes
	for (int j = 0; j < LCM.LinearConstraints();)
	{
		FELinearConstraint& lc = LCM.LinearConstraint(j);

		bool del = false;
		int n = lc.slave.size();
		for (int k = 0; k < n; ++k)
		{
			if (mesh.Node(lc.slave[k].node).HasFlags(FENode::EXCLUDE))
			{
				del = true;
				break;
			}
		}
		if (del) LCM.RemoveLinearConstraint(j); else ++j;
	}

	// reactivate the linear constraints
	LCM.Activate();

	feLog("\tDeactivated elements: %d\n", deactiveElems);
	return (deactiveElems == 0);
}
