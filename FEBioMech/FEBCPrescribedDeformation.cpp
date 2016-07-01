#include "stdafx.h"
#include "FEBCPrescribedDeformation.h"
#include <FECore/FEModel.h>

BEGIN_PARAMETER_LIST(FEBCPrescribedDeformation, FEPrescribedBC)
	ADD_PARAMETER(m_scale, FE_PARAM_DOUBLE  , "scale");
	ADD_PARAMETER(m_F    , FE_PARAM_MAT3D   , "F");
	ADD_PARAMETER(m_G    , FE_PARAM_TENS3DRS, "G");
END_PARAMETER_LIST();

FEBCPrescribedDeformation::FEBCPrescribedDeformation(FEModel* pfem) : FEPrescribedBC(pfem)
{
	m_scale = 1.0;
	m_F.unit();
	m_G.zero();
}

void FEBCPrescribedDeformation::CopyFrom(FEPrescribedBC* pbc)
{
	FEBCPrescribedDeformation* ps = dynamic_cast<FEBCPrescribedDeformation*>(pbc); assert(ps);
	m_scale = ps->m_scale;
	m_F = ps->m_F;
	m_G = ps->m_G;
	m_node = ps->m_node;
	CopyParameterListState(ps->GetParameterList());
}

void FEBCPrescribedDeformation::AddNode(int n)
{
	m_node.push_back(n);
}

void FEBCPrescribedDeformation::AddNodes(const FENodeSet& set)
{
	int N = set.size();
	m_node.resize(N);
	for (int i = 0; i<N; ++i) m_node[i] = set[i];
}

void FEBCPrescribedDeformation::Activate()
{
	FEModel& fem = *GetFEModel();
	DOFS& dofs = fem.GetDOFS();
	int dofX = dofs.GetDOF("x");
	int dofY = dofs.GetDOF("y");
	int dofZ = dofs.GetDOF("z");

	FEMesh& mesh = fem.GetMesh();
	for (size_t i = 0; i<m_node.size(); ++i)
	{
		// get the node
		FENode& node = mesh.Node(m_node[i]);

		// set the dof to prescribed
		node.m_BC[dofX] = DOF_PRESCRIBED;
		node.m_BC[dofY] = DOF_PRESCRIBED;
		node.m_BC[dofZ] = DOF_PRESCRIBED;
	}
}

void FEBCPrescribedDeformation::Deactivate()
{
	FEModel& fem = *GetFEModel();
	DOFS& dofs = fem.GetDOFS();
	int dofX = dofs.GetDOF("x");
	int dofY = dofs.GetDOF("y");
	int dofZ = dofs.GetDOF("z");

	FEMesh& mesh = fem.GetMesh();
	for (size_t i = 0; i<m_node.size(); ++i)
	{
		// get the node
		FENode& node = mesh.Node(m_node[i]);

		// set the dof to prescribed
		node.m_BC[dofX] = DOF_OPEN;
		node.m_BC[dofY] = DOF_OPEN;
		node.m_BC[dofZ] = DOF_OPEN;
	}
}

//-----------------------------------------------------------------------------
void FEBCPrescribedDeformation::SetDeformationGradient(const mat3d& F)
{
	m_F = F;
}

//-----------------------------------------------------------------------------
void FEBCPrescribedDeformation::SetDeformationHessian(const tens3drs& G)
{
	m_G = G;
}

//-----------------------------------------------------------------------------
vec3d FEBCPrescribedDeformation::NodeValue(const vec3d& X)
{
	mat3ds XX = dyad(X);
	vec3d x = m_F*X + m_G.contract2s(XX);
	vec3d u = x - X;
	return u*m_scale;
}

//-----------------------------------------------------------------------------
//! Update the values of the prescribed degrees of freedom.
//! This is called during model update (FESolver::Update)
void FEBCPrescribedDeformation::Update()
{
	FEModel& fem = *GetFEModel();
	DOFS& dofs = fem.GetDOFS();
	int dofX = dofs.GetDOF("x");
	int dofY = dofs.GetDOF("y");
	int dofZ = dofs.GetDOF("z");

	// update the current nodal values
	FEMesh& mesh = fem.GetMesh();
	for (size_t i = 0; i<m_node.size(); ++i)
	{
		FENode& node = mesh.Node(m_node[i]);
		vec3d u = NodeValue(node.m_r0);

		node.set(dofX, u.x);
		node.set(dofY, u.y);
		node.set(dofZ, u.z);
	}
}

//-----------------------------------------------------------------------------
void FEBCPrescribedDeformation::PrepStep(std::vector<double>& ui, bool brel)
{
	FEModel& fem = *GetFEModel();
	DOFS& dofs = fem.GetDOFS();
	int dofX = dofs.GetDOF("x");
	int dofY = dofs.GetDOF("y");
	int dofZ = dofs.GetDOF("z");

	// get the mesh
	FEMesh& mesh = fem.GetMesh();

	int I;
	for (size_t i = 0; i<m_node.size(); ++i)
	{
		FENode& node = mesh.Node(m_node[i]);
		vec3d u = NodeValue(node.m_r0);

		I = -node.m_ID[dofX] - 2; if (I >= 0) ui[I] = (brel ? u.x - node.get(dofX) : u.x);
		I = -node.m_ID[dofY] - 2; if (I >= 0) ui[I] = (brel ? u.y - node.get(dofY) : u.y);
		I = -node.m_ID[dofZ] - 2; if (I >= 0) ui[I] = (brel ? u.z - node.get(dofZ) : u.z);
	}
}
