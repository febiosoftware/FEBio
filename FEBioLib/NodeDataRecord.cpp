#include "stdafx.h"
#include "NodeDataRecord.h"
#include "FECore/FEAnalysis.h"
#include "FEBioMech/FESolidSolver.h"

//-----------------------------------------------------------------------------
double FENodeXPos::value(int nnode) 
{
	FEMesh& mesh = m_pfem->GetMesh();
	FENode& node = mesh.Node(nnode);
	return node.m_rt.x; 
}

//-----------------------------------------------------------------------------
double FENodeYPos::value(int nnode) 
{
	FEMesh& mesh = m_pfem->GetMesh();
	FENode& node = mesh.Node(nnode);
	return node.m_rt.y; 
}

//-----------------------------------------------------------------------------
double FENodeZPos::value(int nnode) 
{
	FEMesh& mesh = m_pfem->GetMesh();
	FENode& node = mesh.Node(nnode);
	return node.m_rt.z; 
}

//-----------------------------------------------------------------------------
double FENodeXDisp::value(int nnode) 
{
	FEMesh& mesh = m_pfem->GetMesh();
	FENode& node = mesh.Node(nnode);
	return node.m_rt.x - node.m_r0.x; 
}

//-----------------------------------------------------------------------------
double FENodeYDisp::value(int nnode) 
{
	FEMesh& mesh = m_pfem->GetMesh();
	FENode& node = mesh.Node(nnode);
	return node.m_rt.y - node.m_r0.y; 
}

//-----------------------------------------------------------------------------
double FENodeZDisp::value(int nnode) 
{
	FEMesh& mesh = m_pfem->GetMesh();
	FENode& node = mesh.Node(nnode);
	return node.m_rt.z - node.m_r0.z; 
}

//-----------------------------------------------------------------------------
double FENodeXVel::value(int nnode) 
{
	FEMesh& mesh = m_pfem->GetMesh();
	FENode& node = mesh.Node(nnode);
	return node.m_vt.x; 
}

//-----------------------------------------------------------------------------
double FENodeYVel::value(int nnode) 
{
	FEMesh& mesh = m_pfem->GetMesh();
	FENode& node = mesh.Node(nnode);
	return node.m_vt.y; 
}

//-----------------------------------------------------------------------------
double FENodeZVel::value(int nnode) 
{
	FEMesh& mesh = m_pfem->GetMesh();
	FENode& node = mesh.Node(nnode);
	return node.m_vt.z; 
}

//-----------------------------------------------------------------------------
double FENodeTemp::value(int nnode) 
{
	FEMesh& mesh = m_pfem->GetMesh();
	FENode& node = mesh.Node(nnode);
	return node.m_T; 
}

//-----------------------------------------------------------------------------
double FENodePressure::value(int nnode) 
{
	FEMesh& mesh = m_pfem->GetMesh();
	FENode& node = mesh.Node(nnode);
	return node.m_pt; 
}

//-----------------------------------------------------------------------------
double FENodeConcentration::value(int nnode) 
{
	FEMesh& mesh = m_pfem->GetMesh();
	FENode& node = mesh.Node(nnode);
	return node.m_ct[0]; 
}

//-----------------------------------------------------------------------------
double FENodeConcentration_::value(int nnode) 
{
	FEMesh& mesh = m_pfem->GetMesh();
	FENode& node = mesh.Node(nnode);
	return node.m_ct[m_nsol]; 
}

//-----------------------------------------------------------------------------
double FENodeForceX::value(int nnode) 
{
	FEMesh& mesh = m_pfem->GetMesh();
	FESolidSolver* psolid_solver = dynamic_cast<FESolidSolver*>(m_pfem->GetCurrentStep()->m_psolver);
	if (psolid_solver)
	{
		vector<double>& Fr = psolid_solver->m_Fr;
		int* id = mesh.Node(nnode).m_ID;
		return (-id[0] - 2 >= 0 ? Fr[-id[0]-2] : 0);
	}
	return 0;
}

//-----------------------------------------------------------------------------
double FENodeForceY::value(int nnode) 
{
	FEMesh& mesh = m_pfem->GetMesh();
	FESolidSolver* psolid_solver = dynamic_cast<FESolidSolver*>(m_pfem->GetCurrentStep()->m_psolver);
	if (psolid_solver)
	{
		vector<double>& Fr = psolid_solver->m_Fr;
		int* id = mesh.Node(nnode).m_ID;
		return (-id[1] - 2 >= 0 ? Fr[-id[1]-2] : 0);
	}
	return 0;
}

//-----------------------------------------------------------------------------
double FENodeForceZ::value(int nnode) 
{
	FEMesh& mesh = m_pfem->GetMesh();
	FESolidSolver* psolid_solver = dynamic_cast<FESolidSolver*>(m_pfem->GetCurrentStep()->m_psolver);
	if (psolid_solver)
	{
		vector<double>& Fr = psolid_solver->m_Fr;
		int* id = mesh.Node(nnode).m_ID;
		return (-id[2] - 2 >= 0 ? Fr[-id[2]-2] : 0);
	}
	return 0;
}

//-----------------------------------------------------------------------------
void NodeDataRecord::Parse(const char* szexpr)
{
	char szcopy[MAX_STRING] = {0};
	strcpy(szcopy, szexpr);
	char* sz = szcopy, *ch;
	m_Data.clear();
	do
	{
		ch = strchr(sz, ';');
		if (ch) *ch++ = 0;
		if      (strcmp(sz, "x" ) == 0) m_Data.push_back(new FENodeXPos(m_pfem));
		else if (strcmp(sz, "y" ) == 0) m_Data.push_back(new FENodeYPos(m_pfem));
		else if (strcmp(sz, "z" ) == 0) m_Data.push_back(new FENodeZPos(m_pfem));
		else if (strcmp(sz, "ux") == 0) m_Data.push_back(new FENodeXDisp(m_pfem));
		else if (strcmp(sz, "uy") == 0) m_Data.push_back(new FENodeYDisp(m_pfem));
		else if (strcmp(sz, "uz") == 0) m_Data.push_back(new FENodeZDisp(m_pfem));
		else if (strcmp(sz, "vx") == 0) m_Data.push_back(new FENodeXVel(m_pfem));
		else if (strcmp(sz, "vy") == 0) m_Data.push_back(new FENodeYVel(m_pfem));
		else if (strcmp(sz, "vz") == 0) m_Data.push_back(new FENodeZVel(m_pfem));
		else if (strcmp(sz, "T" ) == 0) m_Data.push_back(new FENodeTemp(m_pfem));
		else if (strcmp(sz, "p" ) == 0) m_Data.push_back(new FENodePressure(m_pfem));
		else if (strcmp(sz, "c" ) == 0) m_Data.push_back(new FENodeConcentration(m_pfem) );
		else if (strcmp(sz, "c1" ) == 0) m_Data.push_back(new FENodeConcentration_T<0>(m_pfem));
		else if (strcmp(sz, "c2" ) == 0) m_Data.push_back(new FENodeConcentration_T<1>(m_pfem));
		else if (strcmp(sz, "c3" ) == 0) m_Data.push_back(new FENodeConcentration_T<2>(m_pfem));
		else if (strcmp(sz, "c4" ) == 0) m_Data.push_back(new FENodeConcentration_T<3>(m_pfem));
		else if (strcmp(sz, "c5" ) == 0) m_Data.push_back(new FENodeConcentration_T<4>(m_pfem));
		else if (strcmp(sz, "c6" ) == 0) m_Data.push_back(new FENodeConcentration_T<5>(m_pfem));
		else if (strcmp(sz, "c7" ) == 0) m_Data.push_back(new FENodeConcentration_T<6>(m_pfem));
		else if (strcmp(sz, "c8" ) == 0) m_Data.push_back(new FENodeConcentration_T<7>(m_pfem));
		else if (strcmp(sz, "Rx") == 0) m_Data.push_back(new FENodeForceX(m_pfem));
		else if (strcmp(sz, "Ry") == 0) m_Data.push_back(new FENodeForceY(m_pfem));
		else if (strcmp(sz, "Rz") == 0) m_Data.push_back(new FENodeForceZ(m_pfem));
		else throw UnknownDataField(sz);
		sz = ch;
	}
	while (ch);
}

//-----------------------------------------------------------------------------
double NodeDataRecord::Evaluate(int item, int ndata)
{
	FEMesh& mesh = m_pfem->GetMesh();
	int nnode = item - 1;
	if ((nnode < 0) || (nnode >= mesh.Nodes())) return 0;
	return m_Data[ndata]->value(nnode);
}

//-----------------------------------------------------------------------------
void NodeDataRecord::SelectAllItems()
{
	int n = m_pfem->GetMesh().Nodes();
	m_item.resize(n);
	for (int i=0; i<n; ++i) m_item[i] = i+1;
}

//-----------------------------------------------------------------------------
void NodeDataRecord::SetItemList(FENodeSet* pns)
{
	int n = pns->size();
	assert(n);
	m_item.resize(n);
	for (int i=0; i<n; ++i) m_item[i] = (*pns)[i];
}
