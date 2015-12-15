#include "FEBioHeatData.h"
#include "FECore/FEModel.h"

//-----------------------------------------------------------------------------
double FENodeTemp::value(int nnode) 
{
	DOFS& dofs = m_pfem->GetDOFS();
	int dof_T = dofs.GetDOF("t");
	FEMesh& mesh = m_pfem->GetMesh();
	FENode& node = mesh.Node(nnode);
	return node.get(dof_T);
}
