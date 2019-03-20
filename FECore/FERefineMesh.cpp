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
#include "FESurfaceLoad.h"

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

	// update surface loads 
	for (int i = 0; i < fem.SurfaceLoads(); ++i)
	{
		FESurfaceLoad& sl = *fem.SurfaceLoad(i);
		UpdateSurfaceLoad(sl);
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

void FERefineMesh::UpdateSurfaceLoad(FESurfaceLoad& surfLoad)
{
	// look-up table for splitting quads
	const int LUT[4][4] = {
		{0, 4, 8, 7},
		{4, 1, 5, 8},
		{7, 8, 6, 3},
		{8, 5, 2, 6}
	};

	FESurface& surf = surfLoad.GetSurface();

	FEMeshTopo& topo = *m_topo;
	int NF0 = surf.Elements();
	
	// figure out which facets to split
	vector<int> faceList = topo.FaceIndexList(surf);
	assert((int)faceList.size() == NF0);

	// count how many faces to split
	int splitFaces = 0;
	for (int i = 0; i < faceList.size(); ++i)
	{
		int iface = faceList[i];
		if (m_faceList[iface] >= 0) splitFaces++;
	}
	if (splitFaces == 0) return;

	// create a copy of the domain
	FESurface oldSurf(surf);

	// reallocate the domain (Assumes Quad faces!)
	int NF1 = NF0 - splitFaces + 4 * (splitFaces);
	surf.Create(NF1);

	// reinitialize the surface
	int n = 0;
	for (int i = 0; i < NF0; ++i)
	{
		FESurfaceElement& el0 = oldSurf.Element(i);

		int iface = faceList[i];
		if (m_faceList[iface] >= 0)
		{
			const FEFaceList::FACE& face = topo.Face(iface);
			const vector<int>& edge = topo.FaceEdgeList(iface);

			int NL[9];
			NL[0] = face.node[0];
			NL[1] = face.node[1];
			NL[2] = face.node[2];
			NL[3] = face.node[3];
			NL[4] = m_edgeList[edge[0]];
			NL[5] = m_edgeList[edge[1]];
			NL[6] = m_edgeList[edge[2]];
			NL[7] = m_edgeList[edge[3]];
			NL[8] = m_faceList[iface];

			for (int j = 0; j < 4; ++j)
			{
				FESurfaceElement& el1 = surf.Element(n++);
				el1.SetType(FE_QUAD4G4);
				el1.m_node[0] = NL[LUT[j][0]];
				el1.m_node[1] = NL[LUT[j][1]];
				el1.m_node[2] = NL[LUT[j][2]];
				el1.m_node[3] = NL[LUT[j][3]];
			}
		}
		else
		{
			FESurfaceElement& el1 = surf.Element(n++);
			el1 = el0;
		}
	}

	surf.Init();
	surfLoad.SetSurface(&surf);
	if (surfLoad.IsActive()) surfLoad.Activate();
}
