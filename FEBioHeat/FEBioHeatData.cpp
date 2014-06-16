#include "FEBioHeatData.h"
#include "FECore/FEModel.h"

//-----------------------------------------------------------------------------
double FENodeTemp::value(int nnode) 
{
	FEMesh& mesh = m_pfem->GetMesh();
	FENode& node = mesh.Node(nnode);
	return node.m_T; 
}
