#include "stdafx.h"
#include "BC.h"
#include "FEModel.h"
#include "FESolver.h"
#include "FERigidBody.h"
#include "log.h"
#include "DOFS.h"

//-----------------------------------------------------------------------------
FENodalLoad::FENodalLoad(FEModel* pfem) : FEBoundaryCondition(FEBC_ID, pfem)
{
	m_s = 1.0;
	m_lc = -1;
	m_bc = -1;
	m_node = -1;
}

//-----------------------------------------------------------------------------
void FENodalLoad::Serialize(DumpStream& ar)
{
	if (ar.IsShallow()) return;

	FEBoundaryCondition::Serialize(ar);
	if (ar.IsSaving())
	{
		ar << m_bc << m_lc << m_node << m_s;
	}
	else
	{
		ar >> m_bc >> m_lc >> m_node >> m_s;
	}
}

//-----------------------------------------------------------------------------
bool FENodalLoad::Init()
{
	int NLC = GetFEModel()->LoadCurves();
	if ((m_lc < 0)||(m_lc >= NLC))
	{
		felog.printf("ERROR: Invalid loadcurve in nodal load %d\n", GetID());
		return false;
	}
	return true;
}

//-----------------------------------------------------------------------------
//! Return the current value of the nodal load
double FENodalLoad::Value()
{
	FEModel& fem = *GetFEModel();
	return m_s*fem.GetLoadCurve(m_lc)->Value();
}

//-----------------------------------------------------------------------------
FEFixedBC::FEFixedBC(FEModel* pfem) : FEBoundaryCondition(FEBC_ID, pfem)
{ 
	m_dof = -1; 
}

//-----------------------------------------------------------------------------
FEFixedBC::FEFixedBC(FEModel* pfem, int node, int dof) : FEBoundaryCondition(FEBC_ID, pfem)
{ 
	m_node.push_back(node); 
	m_dof = dof; 
}

//-----------------------------------------------------------------------------
void FEFixedBC::AddNode(int node)
{
	m_node.push_back(node);
}

//-----------------------------------------------------------------------------
void FEFixedBC::SetDOF(int dof)
{
	m_dof = dof;
}

//-----------------------------------------------------------------------------
void FEFixedBC::Serialize(DumpStream& ar)
{
	if (ar.IsShallow()) return;

	FEBoundaryCondition::Serialize(ar);
	if (ar.IsSaving())
	{
		ar << m_node << m_dof;
	}
	else
	{
		ar >> m_node >> m_dof;
	}
}

//-----------------------------------------------------------------------------
void FEFixedBC::Activate()
{
	FEBoundaryCondition::Activate();
	if (m_dof >= 0)
	{
		FEMesh& mesh = GetFEModel()->GetMesh();
		int n = (int) m_node.size();
		for (int i=0; i<n; ++i)
		{
			// make sure we only activate open dof's
			vector<int>& BC = mesh.Node(m_node[i]).m_BC;
			if (BC[m_dof] == DOF_OPEN) BC[m_dof] = DOF_FIXED;
		}
	}
}

//-----------------------------------------------------------------------------
void FEFixedBC::Deactivate()
{
	FEBoundaryCondition::Deactivate();
	if (m_dof >= 0)
	{
		FEMesh& mesh = GetFEModel()->GetMesh();
		int n = (int) m_node.size();
		for (int i=0; i<n; ++i)
		{
			vector<int>& BC = mesh.Node(m_node[i]).m_BC;
			BC[m_dof] = DOF_OPEN;
		}
	}
}

//-----------------------------------------------------------------------------
FEPrescribedBC::FEPrescribedBC(FEModel* pfem) : FEBoundaryCondition(FEBC_ID, pfem)
{
	m_lc = -1;
	m_scale = 0.0;
	m_dof = -1;
	m_br = false;
}

//-----------------------------------------------------------------------------
FEPrescribedBC::FEPrescribedBC(FEModel* pfem, const FEPrescribedBC& bc) : FEBoundaryCondition(FEBC_ID, pfem)
{
	m_lc    = bc.m_lc;
	m_scale = bc.m_scale;
	m_dof   = bc.m_dof;
	m_br    = bc.m_br;
	m_item  = bc.m_item;
}

//-----------------------------------------------------------------------------
void FEPrescribedBC::AddNode(int nid, double s)
{
	ITEM item = {nid, s, 0.0};
	m_item.push_back(item);
}

//-----------------------------------------------------------------------------
bool FEPrescribedBC::Init()
{
	// don't forget to call the base class
	if (FEBoundaryCondition::Init() == false) return false;

	// check the load curve ID
	FEModel& fem = *GetFEModel();
	int NLC = fem.LoadCurves();
	if ((m_lc < -1) || (m_lc >= NLC))
	{
		felog.printf("ERROR: Invalid loadcurve in prescribed BC %d\n", GetID());
		return false;
	}

	// make sure this is not a rigid node
	FEMesh& mesh = fem.GetMesh();
	int NN = mesh.Nodes();
	for (size_t i=0; i<m_item.size(); ++i)
	{
		int nid = m_item[i].nid;
		if ((nid < 0) || (nid >= NN)) return false;
		if (mesh.Node(nid).m_rid != -1)
		{
			felog.printf("ERROR: Rigid nodes cannot be prescribed.\n");
			return false;
		}
	}

	return true;
}

//-----------------------------------------------------------------------------
void FEPrescribedBC::Activate()
{
	FEBoundaryCondition::Activate();

	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	for (size_t j = 0; j<m_item.size(); ++j)
	{
		// get the node
		FENode& node = mesh.Node(m_item[j].nid);

		// set the dof to prescribed
		node.m_BC[m_dof] = DOF_PRESCRIBED;

		// evaluate the relative offset
		if (m_br)
		{
			assert(m_dof < node.m_val.size());
			double r = node.get(m_dof);
			m_item[j].ref = r;
		}
	}
}

//-----------------------------------------------------------------------------
void FEPrescribedBC::Deactivate()
{
	FEBoundaryCondition::Deactivate();

	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	for (size_t j = 0; j<m_item.size(); ++j)
	{
		FENode& node = mesh.Node(m_item[j].nid);
		node.m_BC[m_dof] = DOF_OPEN;
	}
}

//-----------------------------------------------------------------------------
double FEPrescribedBC::NodeValue(int n) const
{
	ITEM it = m_item[n];
	double val = m_scale*it.scale;
	if (m_lc >= 0) 
	{
		FEModel& fem = *GetFEModel();
		val *= fem.GetLoadCurve(m_lc)->Value();
	}
	if (m_br) val += it.ref;
	return val;
}

//-----------------------------------------------------------------------------
void FEPrescribedBC::Serialize(DumpStream& ar)
{
	if (ar.IsShallow()) return;

	FEBoundaryCondition::Serialize(ar);
	if (ar.IsSaving())
	{
		ar << m_dof << m_lc << m_item << m_scale << m_br;
	}
	else
	{
		ar >> m_dof >> m_lc >> m_item >> m_scale >> m_br;
	}
}

//-----------------------------------------------------------------------------
//! Update the values of the prescribed degrees of freedom.
//! This is called during model update (FESolver::Update)
void FEPrescribedBC::Update()
{
	// get the mesh
	FEMesh& mesh = GetFEModel()->GetMesh();

	// update the current nodal values
	for (size_t i = 0; i<m_item.size(); ++i)
	{
		FENode& node = mesh.Node(m_item[i].nid);
		double g = NodeValue(i);
		node.set(m_dof, g);
	}
}

//-----------------------------------------------------------------------------
void FEPrescribedBC::PrepStep(std::vector<double>& ui, bool brel)
{
	// get the mesh
	FEMesh& mesh = GetFEModel()->GetMesh();

	for (size_t i = 0; i<m_item.size(); ++i)
	{
		FENode& node = mesh.Node(m_item[i].nid);
		double dq = NodeValue(i);
		int I = -node.m_ID[m_dof] - 2;
		if (I >= 0) ui[I] = (brel ? dq - node.get(m_dof) : dq);
	}
}
