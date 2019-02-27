#include "stdafx.h"
#include "FEPrescribedDOF.h"
#include "FENodeSet.h"
#include "FEModel.h"
#include "fecore_error.h"

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEPrescribedDOF, FEPrescribedBC)
	ADD_PARAMETER(m_scale, "scale");
	ADD_PARAMETER(m_br   , "relative");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEPrescribedDOF::FEPrescribedDOF(FEModel* pfem) : FEPrescribedBC(pfem)
{
	m_scale = 0.0;
	m_dof = -1;
	m_br = false;
}

//-----------------------------------------------------------------------------
FEPrescribedDOF::FEPrescribedDOF(FEModel* pfem, const FEPrescribedDOF& bc) : FEPrescribedBC(pfem)
{
	m_scale = bc.m_scale;
	m_dof = bc.m_dof;
	m_br = bc.m_br;
	m_item = bc.m_item;
	CopyParameterListState(bc.GetParameterList());
}

//-----------------------------------------------------------------------------
// Sets the displacement scale factor. An optional load curve index can be given
// of the load curve that will control the scale factor.
FEPrescribedDOF& FEPrescribedDOF::SetScale(double s, int lc)
{
	m_scale = s;
	if (lc >= 0)
	{
		FEParam* p = FindParameterFromData(&m_scale); assert(p);
		GetFEModel()->AttachLoadController(p, lc);
	}
	return *this;
}

//-----------------------------------------------------------------------------
void FEPrescribedDOF::AddNode(int nid, double s)
{
	ITEM item = { nid, s };
	m_item.push_back(item);
}

//-----------------------------------------------------------------------------
void FEPrescribedDOF::AddNodes(const FENodeSet& nset, double s)
{
	// TODO: Avoid this cast
	FENodeSet* cst = const_cast<FENodeSet*>(&nset);
	m_scale.SetItemList(cst);
	int N = nset.size();
	for (int i = 0; i<N; ++i) AddNode(nset[i], s);
}

//-----------------------------------------------------------------------------
bool FEPrescribedDOF::Init()
{
	// don't forget to call the base class
	if (FEPrescribedBC::Init() == false) return false;

	// make sure this is not a rigid node
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	int NN = mesh.Nodes();
	for (size_t i = 0; i<m_item.size(); ++i)
	{
		int nid = m_item[i].nid;
		if ((nid < 0) || (nid >= NN)) return false;
		if (mesh.Node(nid).m_rid != -1)
		{
			return fecore_error("Rigid nodes cannot be prescribed.");
		}
	}

	return true;
}

//-----------------------------------------------------------------------------
void FEPrescribedDOF::Activate()
{
	FEPrescribedBC::Activate();

	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	for (size_t j = 0; j<m_item.size(); ++j)
	{
		// get the node
		FENode& node = mesh.Node(m_item[j].nid);

		// set the dof to prescribed
		node.set_bc(m_dof, DOF_PRESCRIBED);

		// evaluate the relative offset
		if (m_br)
		{
			assert(m_dof < node.dofs());
			double r = node.get(m_dof);
			m_item[j].ref = r;
		}
	}
}

//-----------------------------------------------------------------------------
void FEPrescribedDOF::Deactivate()
{
	FEPrescribedBC::Deactivate();

	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	for (size_t j = 0; j<m_item.size(); ++j)
	{
		FENode& node = mesh.Node(m_item[j].nid);
		node.set_bc(m_dof, DOF_OPEN);
	}
}

//-----------------------------------------------------------------------------
double FEPrescribedDOF::NodeValue(int n)
{
	FEMesh& mesh = GetFEModel()->GetMesh();
	FENode& node = mesh.Node(n);

	FEMaterialPoint mp;
	mp.m_r0 = node.m_r0;
	mp.m_index = n;

	const ITEM& it = m_item[n];
	double val = m_scale(mp);
	if (m_br) val += it.ref;
	return val;
}

//-----------------------------------------------------------------------------
void FEPrescribedDOF::Serialize(DumpStream& ar)
{
	FEPrescribedBC::Serialize(ar);
	if (ar.IsShallow()) return;
	ar & m_dof & m_item;
}

//-----------------------------------------------------------------------------
void FEPrescribedDOF::CopyFrom(FEPrescribedBC* pbc)
{
	FEPrescribedDOF* ps = dynamic_cast<FEPrescribedDOF*>(pbc); assert(ps);
	m_dof = ps->m_dof;
	m_scale = ps->m_scale;
	m_br = ps->m_br;
	m_item = ps->m_item;
	CopyParameterListState(ps->GetParameterList());
}

//-----------------------------------------------------------------------------
//! Update the values of the prescribed degrees of freedom.
//! This is called during model update (FESolver::Update)
void FEPrescribedDOF::Update()
{
	// get the mesh
	FEMesh& mesh = GetFEModel()->GetMesh();

	// update the current nodal values
	for (size_t i = 0; i<m_item.size(); ++i)
	{
		FENode& node = mesh.Node(m_item[i].nid);
		double g = NodeValue((int)i);
		node.set(m_dof, g);
	}
}

//-----------------------------------------------------------------------------
void FEPrescribedDOF::PrepStep(std::vector<double>& ui, bool brel)
{
	// get the mesh
	FEMesh& mesh = GetFEModel()->GetMesh();

	for (size_t i = 0; i<m_item.size(); ++i)
	{
		FENode& node = mesh.Node(m_item[i].nid);
		double dq = NodeValue((int)i);
		int I = -node.m_ID[m_dof] - 2;
		if (I >= 0) ui[I] = (brel ? dq - node.get(m_dof) : dq);
	}
}
