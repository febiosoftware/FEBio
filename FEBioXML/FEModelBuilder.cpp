#include "stdafx.h"
#include "FEModelBuilder.h"
#include <FECore/FEAnalysis.h>
#include <FECore/BC.h>
#include <FECore/FESurfaceLoad.h>
#include <FECore/FEEdgeLoad.h>
#include <FECore/FEInitialCondition.h>
#include <FECore/FEModelLoad.h>
#include <FECore/FERigidSystem.h>
#include <FECore/RigidBC.h>

//-----------------------------------------------------------------------------
FEModelBuilder::FEModelBuilder(FEModel& fem) : m_fem(fem)
{
	m_pStep = 0;	// zero step pointer
	m_nsteps = 0;	// reset step section counter

	// default element type
	m_ntet4  = FE_TET4G1;
	m_nhex8  = FE_HEX8G8;
	m_ntet10 = FE_TET10G8;
	m_ntet15 = FE_TET15G15;
	m_ntet20 = FE_TET20G15;
	m_ntri6  = FE_TRI6G7;
	m_ntri3  = FE_TRI3G3;
	m_ntri7  = FE_TRI7G7;
	m_ntri10 = FE_TRI10G7;

	// 3-field formulation flags
	m_b3field_hex = true;
	m_b3field_tet = false;

	// UT4 formulation off by default
	m_but4 = false;
}

//-----------------------------------------------------------------------------
void FEModelBuilder::SetModuleName(const std::string& moduleName)
{
	m_moduleName = moduleName;
}

//-----------------------------------------------------------------------------
//! Get the module name
const std::string& FEModelBuilder::GetModuleName() const
{
	return m_moduleName;
}

//-----------------------------------------------------------------------------
FEAnalysis* FEModelBuilder::CreateNewStep()
{
	FEAnalysis* pstep = new FEAnalysis(&m_fem);

	// make sure we have a solver defined
	FESolver* psolver = pstep->GetFESolver();
	if (psolver == 0)
	{
		psolver = BuildSolver(m_fem);
		if (psolver == 0) return 0;
		pstep->SetFESolver(psolver);
	}
	return pstep;
}

//-----------------------------------------------------------------------------
FESolver* FEModelBuilder::BuildSolver(FEModel& fem)
{
	const char* sztype = m_moduleName.c_str();
	FESolver* ps = fecore_new<FESolver>(FESOLVER_ID, sztype, &fem);
	return ps;
}

//-----------------------------------------------------------------------------
void FEModelBuilder::NextStep()
{
	// reset the step pointer
	if (m_nsteps != 0) m_pStep = 0;

	// increase the step section counter
	++m_nsteps;
}

//-----------------------------------------------------------------------------
FEAnalysis* FEModelBuilder::GetStep()
{
	if (m_pStep == 0)
	{
		m_pStep = CreateNewStep();
		m_fem.AddStep(m_pStep);
		if (m_fem.Steps() == 1)
		{
			m_fem.SetCurrentStep(m_pStep);
			m_fem.SetCurrentStepIndex(0);
		}
	}
	return m_pStep;
}

//-----------------------------------------------------------------------------
void FEModelBuilder::AddComponent(FEModelComponent* pmc)
{
	if (m_nsteps > 0)
	{
		GetStep()->AddModelComponent(pmc);
		pmc->Deactivate();
	}
}

//-----------------------------------------------------------------------------
void FEModelBuilder::AddFixedBC(FEFixedBC* pbc)
{
	m_fem.AddFixedBC(pbc);
	AddComponent(pbc);
}

//-----------------------------------------------------------------------------
void FEModelBuilder::AddPrescribedBC(FEPrescribedBC* pbc)
{
	m_fem.AddPrescribedBC(pbc);
	AddComponent(pbc);
}

//-----------------------------------------------------------------------------
void FEModelBuilder::AddNodalLoad(FENodalLoad* pfc)
{
	m_fem.AddNodalLoad(pfc);
	AddComponent(pfc);
}

//-----------------------------------------------------------------------------
void FEModelBuilder::AddSurfaceLoad(FESurfaceLoad* psl)
{
	m_fem.AddSurfaceLoad(psl);
	AddComponent(psl);
}

//-----------------------------------------------------------------------------
void FEModelBuilder::AddEdgeLoad(FEEdgeLoad* pel)
{
	m_fem.AddEdgeLoad(pel);
	AddComponent(pel);
}

//-----------------------------------------------------------------------------
void FEModelBuilder::AddInitialCondition(FEInitialCondition* pic)
{
	m_fem.AddInitialCondition(pic);
	AddComponent(pic);
}

//-----------------------------------------------------------------------------
void FEModelBuilder::AddContactInterface(FESurfacePairInteraction* pci)
{
	m_fem.AddSurfacePairInteraction(pci);
	AddComponent(pci);
}

//-----------------------------------------------------------------------------
void FEModelBuilder::AddModelLoad(FEModelLoad* pml)
{
	m_fem.AddModelLoad(pml);
	AddComponent(pml);
}

//-----------------------------------------------------------------------------
void FEModelBuilder::AddNonlinearConstraint(FENLConstraint* pnc)
{
	m_fem.AddNonlinearConstraint(pnc);
	AddComponent(pnc);
}

//-----------------------------------------------------------------------------
void FEModelBuilder::AddRigidFixedBC(FERigidBodyFixedBC* prc)
{
	m_fem.GetRigidSystem()->AddFixedBC(prc);
	AddComponent(prc);
}

//-----------------------------------------------------------------------------
void FEModelBuilder::AddRigidPrescribedBC(FERigidBodyDisplacement* prc)
{
	m_fem.GetRigidSystem()->AddPrescribedBC(prc);
	AddComponent(prc);	
}

//-----------------------------------------------------------------------------
void FEModelBuilder::AddRigidBodyVelocity(FERigidBodyVelocity* prv)
{
	m_fem.GetRigidSystem()->AddInitialVelocity(prv);
	AddComponent(prv);
}

//-----------------------------------------------------------------------------
void FEModelBuilder::AddRigidBodyAngularVelocity(FERigidBodyAngularVelocity* prv)
{
	m_fem.GetRigidSystem()->AddInitialAngularVelocity(prv);
	AddComponent(prv);
}

//-----------------------------------------------------------------------------
void FEModelBuilder::AddRigidNodeSet(FERigidNodeSet* rs)
{
	m_fem.GetRigidSystem()->AddRigidNodeSet(rs);
	AddComponent(rs);
}

//-----------------------------------------------------------------------------
bool FEModelBuilder::BuildSurface(FESurface& s, FEFacetSet& fs)
{
	FEMesh& m = m_fem.GetMesh();
	int NN = m.Nodes();

	// count nr of faces
	int faces = fs.Faces();

	// allocate storage for faces
	s.Create(faces);

	// read faces
	for (int i = 0; i<faces; ++i)
	{
		FESurfaceElement& el = s.Element(i);
		FEFacetSet::FACET& fi = fs.Face(i);

		if (fi.ntype == 4) el.SetType(FE_QUAD4G4);
		else if (fi.ntype == 3) el.SetType(m_ntri3);
		else if (fi.ntype == 6) el.SetType(m_ntri6);
		else if (fi.ntype == 7) el.SetType(m_ntri7);
		else if (fi.ntype == 10) el.SetType(m_ntri10);
		else if (fi.ntype == 8) el.SetType(FE_QUAD8G9);
		else if (fi.ntype == 9) el.SetType(FE_QUAD9G9);
		else return false;

		int N = el.Nodes(); assert(N == fi.ntype);
		for (int j = 0; j<N; ++j) el.m_node[j] = fi.node[j];
	}

	// copy the name
	s.SetName(fs.GetName());

	return true;
}

//-----------------------------------------------------------------------------
FEModelBuilder::NodeSetPair* FEModelBuilder::FindNodeSetPair(const char* szname)
{
	for (int i = 0; i<m_nsetPair.size(); ++i)
		if (strcmp(m_nsetPair[i].szname, szname) == 0) return &m_nsetPair[i];
	return 0;
}

//-----------------------------------------------------------------------------
FEModelBuilder::NodeSetSet* FEModelBuilder::FindNodeSetSet(const char* szname)
{
	for (int i = 0; i<m_nsetSet.size(); ++i)
		if (strcmp(m_nsetSet[i].szname, szname) == 0) return &m_nsetSet[i];
	return 0;
}

//-----------------------------------------------------------------------------
void FEModelBuilder::BuildNodeList()
{
	// find the min, max ID
	// (We assume that they are given by the first and last node)
	FEMesh& mesh = m_fem.GetMesh();
	int NN = mesh.Nodes();
	int nmin = mesh.Node(0).GetID();
	int nmax = mesh.Node(NN - 1).GetID();
	assert(nmax >= nmin);

	// get the range
	int nn = nmax - nmin + 1;

	// allocate list
	m_node_off = nmin;
	m_node_list.assign(nn, -1);

	// build the list
	for (int i = 0; i<NN; ++i)
	{
		int nid = mesh.Node(i).GetID();
		m_node_list[nid - m_node_off] = i;
	}
}

//-----------------------------------------------------------------------------
int FEModelBuilder::FindNodeFromID(int nid)
{
	int N = (int)m_node_list.size();
	if (N > 0)
	{
		int n = nid - m_node_off;
		if ((n >= 0) && (n<N)) return m_node_list[n];
	}
	return -1;
}

//-----------------------------------------------------------------------------
void FEModelBuilder::GlobalToLocalID(int* l, int n, vector<int>& m)
{
	assert((int)m.size() == n);
	for (int i = 0; i<n; ++i)
	{
		m[i] = FindNodeFromID(l[i]);
	}
}
