#include "stdafx.h"
#include "FEPrescribedNormalDisplacement.h"
#include <FECore/FEModel.h>
#include <FECore/FESurface.h>

BEGIN_PARAMETER_LIST(FEPrescribedNormalDisplacement, FEPrescribedBC)
	ADD_PARAMETER(m_scale, FE_PARAM_DOUBLE, "scale");
	ADD_PARAMETER(m_hint, FE_PARAM_INT, "surface_hint");
END_PARAMETER_LIST()

FEPrescribedNormalDisplacement::FEPrescribedNormalDisplacement(FEModel* fem) : FEPrescribedBC(fem)
{
	m_scale = 0.0;
	m_hint = 0;
}

bool FEPrescribedNormalDisplacement::Init()
{
	return FEPrescribedBC::Init();
}

void FEPrescribedNormalDisplacement::Activate()
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
		FENode& node = mesh.Node(m_node[i].nodeId);

		// set the dof to prescribed
		node.m_BC[dofX] = DOF_PRESCRIBED;
		node.m_BC[dofY] = DOF_PRESCRIBED;
		node.m_BC[dofZ] = DOF_PRESCRIBED;
	}
}

void FEPrescribedNormalDisplacement::Deactivate()
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
		FENode& node = mesh.Node(m_node[i].nodeId);

		// set the dof to prescribed
		node.m_BC[dofX] = DOF_OPEN;
		node.m_BC[dofY] = DOF_OPEN;
		node.m_BC[dofZ] = DOF_OPEN;
	}
}

// assign a node set to the prescribed BC
void FEPrescribedNormalDisplacement::AddNodes(const FESurface& surf)
{
	int N = surf.Nodes();
	m_node.resize(N);
	for (int i=0; i<N; ++i)
	{
		m_node[i].nodeId = surf.NodeIndex(i);
		m_node[i].normal = vec3d(0,0,0);
	}

	if (m_hint == 0)
	{
		int NF = surf.Elements();
		for (int i=0; i<NF; ++i)
		{
			const FESurfaceElement& el = surf.Element(i);
			int nn = el.Nodes();
			if ((nn==3) || (nn == 4))
			{
				for (int n=0; n<nn; ++n)
				{
					int i0 = el.m_lnode[n];
					int ip = el.m_lnode[(n+1)%nn];
					int im = el.m_lnode[(n + nn - 1)%nn];

					vec3d r0 = surf.Node(i0).m_r0;
					vec3d rp = surf.Node(ip).m_r0;
					vec3d rm = surf.Node(im).m_r0;

					vec3d nu = (rp - r0) ^ (rm - r0);

					m_node[i0].normal += nu;
				}
			}
			else if (nn == 6)
			{
				vec3d normals[6];

				// corner nodes
				for (int n = 0; n<3; ++n)
				{
					int i0 = el.m_lnode[n];
					int ip = el.m_lnode[(n + 1) % 3];
					int im = el.m_lnode[(n + 2) % 3];

					vec3d r0 = surf.Node(i0).m_r0;
					vec3d rp = surf.Node(ip).m_r0;
					vec3d rm = surf.Node(im).m_r0;

					vec3d nu = (rp - r0) ^ (rm - r0);

					normals[n] = nu;
				}

				// edge nodes
				for (int n=0; n<3; ++n)
				{
					int n0 = n + 3;
					int n1 = n;
					int n2 = (n+1) % 3;

					normals[n0] = (normals[n1] + normals[n2]) * 0.5;
				}

				for (int n=0; n<6; ++n) m_node[el.m_lnode[n]].normal += normals[n];
			}
		}

		for (int i = 0; i<N; ++i)
		{
			m_node[i].normal.unit();
		}
	}
	else
	{
		int NN = surf.Nodes();
		for (int i=0; i<NN; ++i)
		{
			vec3d ri = -surf.Node(i).m_r0;
			ri.unit();
			m_node[i].normal = ri;
		}
	}
}

// This function is called when the solver needs to know the 
// prescribed dof values. The brel flag indicates wheter the total 
// value is needed or the value with respect to the current nodal dof value
void FEPrescribedNormalDisplacement::PrepStep(std::vector<double>& ui, bool brel)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	DOFS& dofs = fem.GetDOFS();
	int dofX = dofs.GetDOF("x");
	int dofY = dofs.GetDOF("y");
	int dofZ = dofs.GetDOF("z");

	int I;
	for (size_t i = 0; i<m_node.size(); ++i)
	{
		FENode& node = mesh.Node(m_node[i].nodeId);
		vec3d u = m_node[i].normal*m_scale;

		I = -node.m_ID[dofX] - 2; if (I >= 0) ui[I] = u.x;
		I = -node.m_ID[dofY] - 2; if (I >= 0) ui[I] = u.y;
		I = -node.m_ID[dofZ] - 2; if (I >= 0) ui[I] = u.z;
	}
}

// This is called during nodal update and should be used to enforce the 
// nodal degrees of freedoms
void FEPrescribedNormalDisplacement::Update()
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
		FENode& node = mesh.Node(m_node[i].nodeId);
		vec3d u = m_node[i].normal*m_scale;

		node.set(dofX, u.x);
		node.set(dofY, u.y);
		node.set(dofZ, u.z);
	}
}

// copy data from another class
void FEPrescribedNormalDisplacement::CopyFrom(FEPrescribedBC* pbc)
{
	FEPrescribedNormalDisplacement* pnd = dynamic_cast<FEPrescribedNormalDisplacement*>(pbc);
	assert(pnd);
	m_scale = pnd->m_scale;
	m_node = pnd->m_node;
	CopyParameterListState(pnd->GetParameterList());
}
