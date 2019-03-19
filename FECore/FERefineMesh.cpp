#include "stdafx.h"
#include "FERefineMesh.h"
#include "FEModel.h"
#include "FESolidDomain.h"
#include "FEEdgeList.h"
#include "FEElementList.h"
#include "FEFaceList.h"
#include "FEFixedBC.h"
#include "FEPrescribedDOF.h"
#include "FEMeshTopo.h"
#include "FELinearConstraintManager.h"

FERefineMesh::FERefineMesh(FEModel* fem) : FEMeshAdaptor(fem), m_topo(nullptr)
{
}

bool FERefineMesh::BuildMeshTopo()
{
	FEModel& fem = *GetFEModel();
	if (m_topo) { delete m_topo; m_topo = nullptr; }
	m_topo = new FEMeshTopo;
	return m_topo->Create(&fem.GetMesh());
}

void FERefineMesh::UpdateBCs()
{
	FEModel& fem = *GetFEModel();

	// translate BCs
	for (int i = 0; i < fem.FixedBCs(); ++i)
	{
		FEFixedBC& bc = *fem.FixedBC(i);
		UpdateFixedBC(bc);
	}

	for (int i = 0; i < fem.PrescribedBCs(); ++i)
	{
		FEPrescribedDOF& bc = dynamic_cast<FEPrescribedDOF&>(*fem.PrescribedBC(i));
		UpdatePrescribedBC(bc);
	}

	// reactivate the linear constraints
	fem.GetLinearConstraintManager().Activate();
}

void FERefineMesh::UpdateFixedBC(FEFixedBC& bc)
{
	FEMesh& mesh = GetFEModel()->GetMesh();
	vector<int> nodeList = bc.GetNodeList();

	FEMeshTopo& topo = *m_topo;
	vector<int> tag(m_NN, 0);

	for (int j = 0; j < (int)nodeList.size(); ++j) tag[nodeList[j]] = 1;

	for (int j = 0; j < topo.Edges(); ++j)
	{
		if (m_edgeList[j] >= 0)
		{
			const FEEdgeList::EDGE& edge = topo.Edge(j);

			if ((tag[edge.node[0]] == 1) && (tag[edge.node[1]] == 1))
			{
				nodeList.push_back(m_edgeList[j]);
			}
		}
	}

	for (int j = 0; j < topo.Faces(); ++j)
	{
		if (m_faceList[j] >= 0)
		{
			const FEFaceList::FACE& face = topo.Face(j);

			assert(face.ntype == 4);
			if ((tag[face.node[0]] == 1) &&
				(tag[face.node[1]] == 1) &&
				(tag[face.node[2]] == 1) &&
				(tag[face.node[3]] == 1))
			{
				nodeList.push_back(m_faceList[j]);
			}
		}
	}

	// set the node list
	bc.SetNodeList(nodeList);

	// re-activate the bc
	if (bc.IsActive()) bc.Activate();
}

void FERefineMesh::UpdatePrescribedBC(FEPrescribedDOF& bc)
{
	int items = (int)bc.Items();

	FEMeshTopo& topo = *m_topo;
	vector<int> tag(m_NN, -1);

	for (int j = 0; j < items; ++j) tag[bc.GetItem(j).nid] = j;

	for (int j = 0; j < topo.Edges(); ++j)
	{
		if (m_edgeList[j] >= 0)
		{
			const FEEdgeList::EDGE& edge = topo.Edge(j);

			if ((tag[edge.node[0]] >= 0) && (tag[edge.node[1]] >= 0))
			{
				double a0 = bc.GetItem(tag[edge.node[0]]).ref;
				double a1 = bc.GetItem(tag[edge.node[1]]).ref;
				bc.AddNode(m_edgeList[j], (a0 + a1)*0.5);
			}
		}
	}

	for (int j = 0; j < topo.Faces(); ++j)
	{
		if (m_faceList[j] >= 0)
		{
			const FEFaceList::FACE& face = topo.Face(j);

			assert(face.ntype == 4);
			if ((tag[face.node[0]] >= 0) &&
				(tag[face.node[1]] >= 0) &&
				(tag[face.node[2]] >= 0) &&
				(tag[face.node[3]] >= 0))
			{
				double a0 = bc.GetItem(tag[face.node[0]]).ref;
				double a1 = bc.GetItem(tag[face.node[1]]).ref;
				double a2 = bc.GetItem(tag[face.node[2]]).ref;
				double a3 = bc.GetItem(tag[face.node[3]]).ref;
				bc.AddNode(m_faceList[j], (a0 + a1 + a2 + a3)*0.25);
			}
		}
	}

	// re-activate the bc
	if (bc.IsActive()) bc.Activate();
}
