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

FERefineMesh::FERefineMesh(FEModel* fem) : FEMeshAdaptor(fem), m_topo(nullptr)
{
}

bool FERefineMesh::BuildMeshTopo(FEModel& fem)
{
	if (m_topo) { delete m_topo; m_topo = nullptr; }
	m_topo = new FEMeshTopo;
	return m_topo->Create(&fem.GetMesh());
}

void FERefineMesh::UpdateFixedBC(FEFixedBC& bc)
{
	FEMesh& mesh = GetFEModel()->GetMesh();
	vector<int> nodeList = bc.GetNodeList();

	FEMeshTopo& topo = *m_topo;
	vector<int> tag(m_NN, 0);

	for (int j = 0; j < (int)nodeList.size(); ++j) tag[nodeList[j]] = 1;

	for (int j = 0; j < topo.m_edgeList.Edges(); ++j)
	{
		const FEEdgeList::EDGE& edge = topo.m_edgeList[j];

		if ((tag[edge.node[0]] == 1) && (tag[edge.node[1]] == 1))
		{
			nodeList.push_back(m_N0 + j);
		}
	}

	for (int j = 0; j < topo.m_faceList.Faces(); ++j)
	{
		const FEFaceList::FACE& face = topo.m_faceList.Face(j);

		assert(face.ntype == 4);
		if ((tag[face.node[0]] == 1) &&
			(tag[face.node[1]] == 1) &&
			(tag[face.node[2]] == 1) &&
			(tag[face.node[3]] == 1))
		{
			nodeList.push_back(m_N0 + m_NC + j);
		}
	}

	// set the node list
	bc.SetNodeList(nodeList);

	// re-activate the bc
	if (bc.IsActive()) bc.Activate();
}

void FERefineMesh::UpdatePrescribedBC(FEPrescribedDOF& bc)
{
	int items = bc.Items();

	FEMeshTopo& topo = *m_topo;
	vector<int> tag(m_NN, -1);

	for (int j = 0; j < items; ++j) tag[bc.GetItem(j).nid] = j;

	for (int j = 0; j < topo.m_edgeList.Edges(); ++j)
	{
		const FEEdgeList::EDGE& edge = topo.m_edgeList[j];

		if ((tag[edge.node[0]] >= 0) && (tag[edge.node[1]] >= 0))
		{
			double a0 = bc.GetItem(tag[edge.node[0]]).ref;
			double a1 = bc.GetItem(tag[edge.node[1]]).ref;
			bc.AddNode(m_N0 + j, (a0 + a1)*0.5);
		}
	}

	for (int j = 0; j < topo.m_faceList.Faces(); ++j)
	{
		const FEFaceList::FACE& face = topo.m_faceList.Face(j);

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
			bc.AddNode(m_N0 + m_NC + j, (a0 + a1 + a2 + a3)*0.25);
		}
	}

	// re-activate the bc
	if (bc.IsActive()) bc.Activate();
}
