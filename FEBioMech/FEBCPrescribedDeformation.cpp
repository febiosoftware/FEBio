#include "stdafx.h"
#include "FEBCPrescribedDeformation.h"
#include <FECore/FELinearConstraintManager.h>
#include <FECore/FEModel.h>

BEGIN_PARAMETER_LIST(FEBCPrescribedDeformation, FEPrescribedBC)
	ADD_PARAMETER(m_scale, "scale");
	ADD_PARAMETER(m_F    , "F");
END_PARAMETER_LIST();

FEBCPrescribedDeformation::FEBCPrescribedDeformation(FEModel* pfem) : FEPrescribedBC(pfem)
{
	m_scale = 1.0;
	m_F.unit();
}

void FEBCPrescribedDeformation::CopyFrom(FEPrescribedBC* pbc)
{
	FEBCPrescribedDeformation* ps = dynamic_cast<FEBCPrescribedDeformation*>(pbc); assert(ps);
	m_scale = ps->m_scale;
	m_F = ps->m_F;
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
vec3d FEBCPrescribedDeformation::NodeValue(const vec3d& X)
{
	mat3ds XX = dyad(X);
	vec3d x = m_F*X;
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

//=============================================================================
BEGIN_PARAMETER_LIST(FEBCPrescribedDeformation2O, FEPrescribedBC)
	ADD_PARAMETER(m_scale, "scale");
	ADD_PARAMETER(m_F    , "F");
	ADD_PARAMETER(m_G    , "G");
	ADD_PARAMETER(m_refNode, "reference");
END_PARAMETER_LIST();

FEBCPrescribedDeformation2O::FEBCPrescribedDeformation2O(FEModel* pfem) : FEPrescribedBC(pfem)
{
	m_scale = 1.0;
	m_F.unit();
	m_G.zero();
	m_refNode = -1;
}

bool FEBCPrescribedDeformation2O::Init()
{
	if (FEPrescribedBC::Init() == false) return false;

//	m_refNode--;
	if (m_refNode < 0) return false;

	return true;
}

//-----------------------------------------------------------------------------
// Sets the displacement scale factor. An optional load curve index can be given
// of the load curve that will control the scale factor.
void FEBCPrescribedDeformation2O::SetScale(double s, int lc)
{
	m_scale = s;
	if (lc >= 0)
	{
		FEParam& p = *FEParamContainer::FindParameterFromData((void*)(&m_scale));
		p.SetLoadCurve(lc, m_scale);
	}
}

void FEBCPrescribedDeformation2O::CopyFrom(FEPrescribedBC* pbc)
{
	FEBCPrescribedDeformation2O* ps = dynamic_cast<FEBCPrescribedDeformation2O*>(pbc); assert(ps);
	m_scale = ps->m_scale;
	m_F = ps->m_F;
	m_G = ps->m_G;
	m_node = ps->m_node;
	m_refNode = ps->m_refNode;
	CopyParameterListState(ps->GetParameterList());
}

void FEBCPrescribedDeformation2O::SetReferenceNode(int n)
{
	m_refNode = n;
}

void FEBCPrescribedDeformation2O::AddNode(int n)
{
	m_node.push_back(n);
}

void FEBCPrescribedDeformation2O::AddNodes(const FENodeSet& set)
{
	int N = set.size();
	m_node.resize(N);
	for (int i = 0; i<N; ++i) m_node[i] = set[i];
}

void FEBCPrescribedDeformation2O::Activate()
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

void FEBCPrescribedDeformation2O::Deactivate()
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
void FEBCPrescribedDeformation2O::SetDeformationGradient(const mat3d& F)
{
	m_F = F;
}

//-----------------------------------------------------------------------------
void FEBCPrescribedDeformation2O::SetDeformationHessian(const tens3drs& G)
{
	m_G = G;
}

//-----------------------------------------------------------------------------
vec3d FEBCPrescribedDeformation2O::NodeValue(const vec3d& X1, const vec3d& X)
{
	mat3ds XX = dyad(X);
	mat3ds XX1 = dyad(X1);
	mat3d U = m_F - mat3dd(1.0);
	vec3d u = U*(X - X1) + m_G.contract2s(XX - XX1)*0.5;
	return u*m_scale;
}

//-----------------------------------------------------------------------------
//! Update the values of the prescribed degrees of freedom.
//! This is called during model update (FESolver::Update)
void FEBCPrescribedDeformation2O::Update()
{
	FEModel& fem = *GetFEModel();
	DOFS& dofs = fem.GetDOFS();
	int dofX = dofs.GetDOF("x");
	int dofY = dofs.GetDOF("y");
	int dofZ = dofs.GetDOF("z");

	// update the current nodal values
	FEMesh& mesh = fem.GetMesh();
	vec3d X1 = mesh.Node(m_refNode).m_r0;
	for (size_t i = 0; i<m_node.size(); ++i)
	{
		FENode& node = mesh.Node(m_node[i]);
		vec3d u = NodeValue(X1, node.m_r0);

		node.set(dofX, u.x);
		node.set(dofY, u.y);
		node.set(dofZ, u.z);
	}
}

//-----------------------------------------------------------------------------
void FEBCPrescribedDeformation2O::PrepStep(std::vector<double>& ui, bool brel)
{
	FEModel& fem = *GetFEModel();
	DOFS& dofs = fem.GetDOFS();
	int dofX = dofs.GetDOF("x");
	int dofY = dofs.GetDOF("y");
	int dofZ = dofs.GetDOF("z");

	// get the mesh
	FEMesh& mesh = fem.GetMesh();
	vec3d X1 = mesh.Node(m_refNode).m_r0;

	int I;
	for (size_t i = 0; i<m_node.size(); ++i)
	{
		FENode& node = mesh.Node(m_node[i]);
		vec3d u = NodeValue(X1, node.m_r0);

		I = -node.m_ID[dofX] - 2; if (I >= 0) ui[I] = (brel ? u.x - node.get(dofX) : u.x);
		I = -node.m_ID[dofY] - 2; if (I >= 0) ui[I] = (brel ? u.y - node.get(dofY) : u.y);
		I = -node.m_ID[dofZ] - 2; if (I >= 0) ui[I] = (brel ? u.z - node.get(dofZ) : u.z);
	}


	// THIS IS A HACK!!!! FIX THIS!!!
/*
	// get the mesh
	mat3d U = m_F - mat3dd(1.0);
	FELinearConstraintManager& LM = fem.GetLinearConstraintManager();
	const int NL = LM.LinearConstraints();
	for (int i = 0; i<NL; ++i)
	{
		FELinearConstraint& lc = LM.LinearConstraint(i);

		FENode& slaveNode = mesh.Node(lc.master.node);
		FENode& masterNode = mesh.Node(lc.slave[0].node);

		vec3d Xp = slaveNode.m_r0;
		vec3d Xm = masterNode.m_r0;

		mat3ds XXp = dyad(Xp);
		mat3ds XXm = dyad(Xm);

		vec3d d = (U*(Xp - Xm) + (m_G.contract2s(XXp - XXm))*0.5)*m_scale;

		switch (lc.master.dof)
		{
		case 0: lc.m_off = d.x; break;
		case 1: lc.m_off = d.y; break;
		case 2: lc.m_off = d.z; break;
		}

	}
*/
	// END HACK
}
