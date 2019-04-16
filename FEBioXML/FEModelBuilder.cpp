/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, Columbia University, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/

#include "stdafx.h"
#include "FEModelBuilder.h"
#include <FECore/FEMaterial.h>
#include <FECore/FEAnalysis.h>
#include <FECore/FEBoundaryCondition.h>
#include <FECore/FENodalLoad.h>
#include <FECore/FESurfaceLoad.h>
#include <FECore/FEEdgeLoad.h>
#include <FECore/FEInitialCondition.h>
#include <FECore/FEModelLoad.h>
#include <FEBioMech/FERigidSystem.h>
#include <FEBioMech/RigidBC.h>
#include <FEBioMech/FEUDGHexDomain.h>
#include <FEBioMech/FEUT4Domain.h>
#include <FEBioMech/FEMechModel.h>
#include <FECore/FEEdge.h>

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
    m_b3field_shell = false;
    m_b3field_quad = true;
    m_b3field_tri = false;

	// shell formulation
	m_default_shell = NEW_SHELL;

	// UT4 formulation off by default
	m_but4 = false;
	m_ut4_alpha = 0.05;
	m_ut4_bdev = false;

	// UDG hourglass parameter
	m_udghex_hg = 1.0;
}

//-----------------------------------------------------------------------------
void FEModelBuilder::SetModuleName(const std::string& moduleName)
{
	m_fem.SetModuleName(moduleName);
	FECoreKernel::GetInstance().SetActiveModule(moduleName.c_str());
}

//-----------------------------------------------------------------------------
//! Get the module name
std::string FEModelBuilder::GetModuleName() const
{
	return m_fem.GetModuleName();
}

//-----------------------------------------------------------------------------
FEAnalysis* FEModelBuilder::CreateNewStep()
{
	FEAnalysis* pstep = fecore_new<FEAnalysis>("analysis", &m_fem);

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
// create a material
FEMaterial* FEModelBuilder::CreateMaterial(const char* sztype)
{
	FEMaterial* pmat = fecore_new<FEMaterial>(sztype, &m_fem);
	return pmat;
}

//-----------------------------------------------------------------------------
FESolver* FEModelBuilder::BuildSolver(FEModel& fem)
{
	string moduleName = fem.GetModuleName();
	const char* sztype = moduleName.c_str();
	FESolver* ps = fecore_new<FESolver>(sztype, &fem);
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
//! Create a domain
FEDomain* FEModelBuilder::CreateDomain(FE_Element_Spec espec, FEMaterial* mat)
{
	FECoreKernel& febio = FECoreKernel::GetInstance();
	FEDomain* pdom = febio.CreateDomain(espec, &m_fem.GetMesh(), mat);

	// Handle dome special cases
	// TODO: Find a better way of dealing with these special cases
	FEUDGHexDomain* udg = dynamic_cast<FEUDGHexDomain*>(pdom);
	if (udg)
	{
		udg->SetHourGlassParameter(m_udghex_hg);
	}

	FEUT4Domain* ut4 = dynamic_cast<FEUT4Domain*>(pdom);
	if (ut4)
	{
		ut4->SetUT4Parameters(m_ut4_alpha, m_ut4_bdev);
	}

	return pdom;
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
void FEModelBuilder::AddBC(FEBoundaryCondition* pbc)
{
	m_fem.AddBoundaryCondition(pbc);
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
void FEModelBuilder::AddContactInterface(FESurfacePairConstraint* pci)
{
	m_fem.AddSurfacePairConstraint(pci);
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
	static_cast<FEMechModel&>(m_fem).GetRigidSystem()->AddFixedBC(prc);
	AddComponent(prc);
}

//-----------------------------------------------------------------------------
void FEModelBuilder::AddRigidPrescribedBC(FERigidBodyDisplacement* prc)
{
	static_cast<FEMechModel&>(m_fem).GetRigidSystem()->AddPrescribedBC(prc);
	AddComponent(prc);	
}

//-----------------------------------------------------------------------------
void FEModelBuilder::AddRigidBodyVelocity(FERigidBodyVelocity* prv)
{
	static_cast<FEMechModel&>(m_fem).GetRigidSystem()->AddInitialVelocity(prv);
	AddComponent(prv);
}

//-----------------------------------------------------------------------------
void FEModelBuilder::AddRigidBodyAngularVelocity(FERigidBodyAngularVelocity* prv)
{
	static_cast<FEMechModel&>(m_fem).GetRigidSystem()->AddInitialAngularVelocity(prv);
	AddComponent(prv);
}

//-----------------------------------------------------------------------------
void FEModelBuilder::AddRigidNodeSet(FERigidNodeSet* rs)
{
	static_cast<FEMechModel&>(m_fem).GetRigidSystem()->AddRigidNodeSet(rs);
	AddComponent(rs);
}

//---------------------------------------------------------------------------------
// parse a surface section for contact definitions
//
bool FEModelBuilder::BuildSurface(FESurface& s, FEFacetSet& fs, bool bnodal)
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

		// set the element type/integration rule
		if (bnodal)
		{
			if (fi.ntype == 4) el.SetType(FE_QUAD4NI);
			else if (fi.ntype == 3) el.SetType(FE_TRI3NI);
			else if (fi.ntype == 6) el.SetType(FE_TRI6NI);
			else if (fi.ntype == 8) el.SetType(FE_QUAD8NI);
			else if (fi.ntype == 9) el.SetType(FE_QUAD9NI);
			else return false;
		}
		else
		{
			if (fi.ntype == 4) el.SetType(FE_QUAD4G4);
			else if (fi.ntype == 3) el.SetType(m_ntri3);
			else if (fi.ntype == 6) el.SetType(m_ntri6);
			else if (fi.ntype == 7) el.SetType(m_ntri7);
			else if (fi.ntype == 10) el.SetType(m_ntri10);
			else if (fi.ntype == 8) el.SetType(FE_QUAD8G9);
			else if (fi.ntype == 9) el.SetType(FE_QUAD9G9);
			else return false;
		}

		int N = el.Nodes(); assert(N == fi.ntype);
		for (int j = 0; j<N; ++j) el.m_node[j] = fi.node[j];
	}

    // copy the name
    s.SetName(fs.GetName());

	s.Init();
    
	return true;
}

//-----------------------------------------------------------------------------
bool FEModelBuilder::BuildEdge(FEEdge& e, FESegmentSet& es)
{
	FEMesh& m = m_fem.GetMesh();
	int NN = m.Nodes();

	// count nr of segments
	int nsegs = es.Segments();

	// allocate storage for faces
	e.Create(nsegs);

	// read segments
	for (int i = 0; i<nsegs; ++i)
	{
		FELineElement& el = e.Element(i);
		FESegmentSet::SEGMENT& si = es.Segment(i);

		if (si.ntype == 2) el.SetType(FE_LINE2G1);
		else return false;

		int N = el.Nodes(); assert(N == si.ntype);
		for (int j = 0; j<N; ++j) el.m_node[j] = si.node[j];
	}

	// copy the name
	e.SetName(es.GetName());

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

//-----------------------------------------------------------------------------
// Call this to initialize default variables when reading older files.
void FEModelBuilder::SetDefaultVariables()
{
	// Reset degrees of
	FEModel& fem = m_fem;
	DOFS& dofs = fem.GetDOFS();
	dofs.Reset();

	// Add the default variables and degrees of freedom
	int varD = dofs.AddVariable("displacement", VAR_VEC3);
	dofs.SetDOFName(varD, 0, "x");
	dofs.SetDOFName(varD, 1, "y");
	dofs.SetDOFName(varD, 2, "z");
	int varQ = dofs.AddVariable("rotation", VAR_VEC3);
	dofs.SetDOFName(varQ, 0, "u");
	dofs.SetDOFName(varQ, 1, "v");
	dofs.SetDOFName(varQ, 2, "w");
	int varSD = dofs.AddVariable("shell displacement", VAR_VEC3);
	dofs.SetDOFName(varSD, 0, "sx");
	dofs.SetDOFName(varSD, 1, "sy");
	dofs.SetDOFName(varSD, 2, "sz");
	int varP = dofs.AddVariable("fluid pressure");
	dofs.SetDOFName(varP, 0, "p");
	int varSP = dofs.AddVariable("shell fluid pressure");
	dofs.SetDOFName(varSP, 0, "q");
	int varQR = dofs.AddVariable("rigid rotation", VAR_VEC3);
	dofs.SetDOFName(varQR, 0, "Ru");
	dofs.SetDOFName(varQR, 1, "Rv");
	dofs.SetDOFName(varQR, 2, "Rw");
	int varT = dofs.AddVariable("temperature");
	dofs.SetDOFName(varT, 0, "T");
	int varV = dofs.AddVariable("velocity", VAR_VEC3);
	dofs.SetDOFName(varV, 0, "vx");
	dofs.SetDOFName(varV, 1, "vy");
	dofs.SetDOFName(varV, 2, "vz");
	int varW = dofs.AddVariable("relative fluid velocity", VAR_VEC3);
	dofs.SetDOFName(varW, 0, "wx");
	dofs.SetDOFName(varW, 1, "wy");
	dofs.SetDOFName(varW, 2, "wz");
	int varWP = dofs.AddVariable("previous relative fluid velocity", VAR_VEC3);
	dofs.SetDOFName(varWP, 0, "wxp");
	dofs.SetDOFName(varWP, 1, "wyp");
	dofs.SetDOFName(varWP, 2, "wzp");
	int varAW = dofs.AddVariable("relative fluid acceleration", VAR_VEC3);
	dofs.SetDOFName(varAW, 0, "awx");
	dofs.SetDOFName(varAW, 1, "awy");
	dofs.SetDOFName(varAW, 2, "awz");
	int varAWP = dofs.AddVariable("previous relative fluid acceleration", VAR_VEC3);
	dofs.SetDOFName(varAWP, 0, "awxp");
	dofs.SetDOFName(varAWP, 1, "awyp");
	dofs.SetDOFName(varAWP, 2, "awzp");
	int varVF = dofs.AddVariable("fluid velocity", VAR_VEC3);
	dofs.SetDOFName(varVF, 0, "vfx");
	dofs.SetDOFName(varVF, 1, "vfy");
	dofs.SetDOFName(varVF, 2, "vfz");
	int varAF = dofs.AddVariable("fluid acceleration", VAR_VEC3);
	dofs.SetDOFName(varAF, 0, "afx");
	dofs.SetDOFName(varAF, 1, "afy");
	dofs.SetDOFName(varAF, 2, "afz");
	int varEF = dofs.AddVariable("fluid dilation");
	dofs.SetDOFName(varEF, 0, "ef");
	int varEFP = dofs.AddVariable("previous fluid dilation");
	dofs.SetDOFName(varEFP, 0, "efp");
	int varAEF = dofs.AddVariable("fluid dilation tderiv");
	dofs.SetDOFName(varAEF, 0, "aef");
	int varAEP = dofs.AddVariable("previous fluid dilation tderiv");
	dofs.SetDOFName(varAEP, 0, "aefp");
	int varQP = dofs.AddVariable("previous rotation", VAR_VEC3);
	dofs.SetDOFName(varQP, 0, "up");
	dofs.SetDOFName(varQP, 1, "vp");
	dofs.SetDOFName(varQP, 2, "wp");
	int varSDP = dofs.AddVariable("previous shell displacement", VAR_VEC3);
	dofs.SetDOFName(varSDP, 0, "sxp");
	dofs.SetDOFName(varSDP, 1, "syp");
	dofs.SetDOFName(varSDP, 2, "szp");
	int varQV = dofs.AddVariable("shell velocity", VAR_VEC3);
	dofs.SetDOFName(varQV, 0, "svx");
	dofs.SetDOFName(varQV, 1, "svy");
	dofs.SetDOFName(varQV, 2, "svz");
	int varQA = dofs.AddVariable("shell acceleration", VAR_VEC3);
	dofs.SetDOFName(varQA, 0, "sax");
	dofs.SetDOFName(varQA, 1, "say");
	dofs.SetDOFName(varQA, 2, "saz");
	int varQVP = dofs.AddVariable("previous shell velocity", VAR_VEC3);
	dofs.SetDOFName(varQVP, 0, "svxp");
	dofs.SetDOFName(varQVP, 1, "svyp");
	dofs.SetDOFName(varQVP, 2, "svzp");
	int varQAP = dofs.AddVariable("previous shell acceleration", VAR_VEC3);
	dofs.SetDOFName(varQAP, 0, "saxp");
	dofs.SetDOFName(varQAP, 1, "sayp");
	dofs.SetDOFName(varQAP, 2, "sazp");
	// must be last variable definition!!
	int varC = dofs.AddVariable("concentration", VAR_ARRAY); // we start with zero concentrations
															 // must be last variable definition!!
	int varSC = dofs.AddVariable("shell concentration", VAR_ARRAY); // we start with zero concentrations
}

//-----------------------------------------------------------------------------
//! Get the element type from a XML tag
FE_Element_Spec FEModelBuilder::ElementSpec(const char* sztype)
{
    FEMesh& mesh = m_fem.GetMesh();
    
	// determine the element shape
	FE_Element_Shape eshape = FE_ELEM_INVALID_SHAPE;

	// for shells, don't overwrite m_pim->m_ntri3/6 or m_nquad4/8, since they are needed for surface definitions
	FE_Element_Type stype = FE_ELEM_INVALID_TYPE;
	if      (strcmp(sztype, "hex8"   ) == 0) eshape = ET_HEX8;
	else if (strcmp(sztype, "hex20"  ) == 0) eshape = ET_HEX20;
	else if (strcmp(sztype, "hex27"  ) == 0) eshape = ET_HEX27;
	else if (strcmp(sztype, "penta6" ) == 0) eshape = ET_PENTA6;
	else if (strcmp(sztype, "penta15") == 0) eshape = ET_PENTA15;
	else if (strcmp(sztype, "pyra5"  ) == 0) eshape = ET_PYRA5;
	else if (strcmp(sztype, "tet4"   ) == 0) eshape = ET_TET4;
	else if (strcmp(sztype, "tet10"  ) == 0) eshape = ET_TET10;
	else if (strcmp(sztype, "tet15"  ) == 0) eshape = ET_TET15;
	else if (strcmp(sztype, "tet20"  ) == 0) eshape = ET_TET20;
	else if (strcmp(sztype, "quad4"  ) == 0) { eshape = ET_QUAD4; stype = FE_SHELL_QUAD4G8; }   // default shell type for quad4
	else if (strcmp(sztype, "quad8"  ) == 0) { eshape = ET_QUAD8; stype = FE_SHELL_QUAD8G18; }   // default shell type for quad8
	else if (strcmp(sztype, "quad9"  ) == 0) eshape = ET_QUAD9;
	else if (strcmp(sztype, "tri3"   ) == 0) { eshape = ET_TRI3; stype = FE_SHELL_TRI3G6; }     // default shell type for tri3
	else if (strcmp(sztype, "tri6"   ) == 0) { eshape = ET_TRI6; stype = FE_SHELL_TRI6G14; }     // default shell type for tri6
    else if (strcmp(sztype, "q4eas"  ) == 0) { eshape = ET_QUAD4; stype = FE_SHELL_QUAD4G8; m_default_shell = EAS_SHELL; }   // default shell type for q4eas
    else if (strcmp(sztype, "q4ans"  ) == 0) { eshape = ET_QUAD4; stype = FE_SHELL_QUAD4G8; m_default_shell = ANS_SHELL; }   // default shell type for q4ans
	else if (strcmp(sztype, "truss2" ) == 0) eshape = ET_TRUSS2;
	else
	{
		// new way for defining element type and integration rule at the same time
		// this is useful for multi-step analyses where the geometry is read in before the control section.
		if      (strcmp(sztype, "TET10G4"     ) == 0) { eshape = ET_TET10; m_ntet10 = FE_TET10G4; }
		else if (strcmp(sztype, "TET10G8"     ) == 0) { eshape = ET_TET10; m_ntet10 = FE_TET10G8; }
		else if (strcmp(sztype, "TET10GL11"   ) == 0) { eshape = ET_TET10; m_ntet10 = FE_TET10GL11; }
		else if (strcmp(sztype, "TET10G4_S3"  ) == 0) { eshape = ET_TET10; m_ntet10 = FE_TET10G4;   m_ntri6 = FE_TRI6G3; }
		else if (strcmp(sztype, "TET10G8_S3"  ) == 0) { eshape = ET_TET10; m_ntet10 = FE_TET10G8;   m_ntri6 = FE_TRI6G3; }
		else if (strcmp(sztype, "TET10GL11_S3") == 0) { eshape = ET_TET10; m_ntet10 = FE_TET10GL11; m_ntri6 = FE_TRI6G3; }
		else if (strcmp(sztype, "TET10G4_S4"  ) == 0) { eshape = ET_TET10; m_ntet10 = FE_TET10G4;   m_ntri6 = FE_TRI6G4; }
		else if (strcmp(sztype, "TET10G8_S4"  ) == 0) { eshape = ET_TET10; m_ntet10 = FE_TET10G8;   m_ntri6 = FE_TRI6G4; }
		else if (strcmp(sztype, "TET10GL11_S4") == 0) { eshape = ET_TET10; m_ntet10 = FE_TET10GL11; m_ntri6 = FE_TRI6G4; }
		else if (strcmp(sztype, "TET10G4_S7"  ) == 0) { eshape = ET_TET10; m_ntet10 = FE_TET10G4;   m_ntri6 = FE_TRI6G7; }
		else if (strcmp(sztype, "TET10G8_S7"  ) == 0) { eshape = ET_TET10; m_ntet10 = FE_TET10G8;   m_ntri6 = FE_TRI6G7; }
		else if (strcmp(sztype, "TET10GL11_S7") == 0) { eshape = ET_TET10; m_ntet10 = FE_TET10GL11; m_ntri6 = FE_TRI6G7; }
		else if (strcmp(sztype, "TET15G8"     ) == 0) { eshape = ET_TET15; m_ntet15 = FE_TET15G8; }
		else if (strcmp(sztype, "TET15G11"    ) == 0) { eshape = ET_TET15; m_ntet15 = FE_TET15G11; }
		else if (strcmp(sztype, "TET15G15"    ) == 0) { eshape = ET_TET15; m_ntet15 = FE_TET15G15; }
		else if (strcmp(sztype, "TET15G8_S3"  ) == 0) { eshape = ET_TET15; m_ntet15 = FE_TET15G8;  m_ntri7 = FE_TRI7G3; }
		else if (strcmp(sztype, "TET15G11_S3" ) == 0) { eshape = ET_TET15; m_ntet15 = FE_TET15G11; m_ntri7 = FE_TRI7G3; }
		else if (strcmp(sztype, "TET15G15_S3" ) == 0) { eshape = ET_TET15; m_ntet15 = FE_TET15G15; m_ntri7 = FE_TRI7G3; }
		else if (strcmp(sztype, "TET15G8_S4"  ) == 0) { eshape = ET_TET15; m_ntet15 = FE_TET15G8;  m_ntri7 = FE_TRI7G4; }
		else if (strcmp(sztype, "TET15G11_S4" ) == 0) { eshape = ET_TET15; m_ntet15 = FE_TET15G11; m_ntri7 = FE_TRI7G4; }
		else if (strcmp(sztype, "TET15G15_S4" ) == 0) { eshape = ET_TET15; m_ntet15 = FE_TET15G15; m_ntri7 = FE_TRI7G4; }
		else if (strcmp(sztype, "TET15G8_S7"  ) == 0) { eshape = ET_TET15; m_ntet15 = FE_TET15G8;  m_ntri7 = FE_TRI7G7; }
		else if (strcmp(sztype, "TET15G11_S7" ) == 0) { eshape = ET_TET15; m_ntet15 = FE_TET15G11; m_ntri7 = FE_TRI7G7; }
		else if (strcmp(sztype, "TET15G15_S7" ) == 0) { eshape = ET_TET15; m_ntet15 = FE_TET15G15; m_ntri7 = FE_TRI7G7; }
		else if (strcmp(sztype, "PENTA15G8"   ) == 0) { eshape = ET_PENTA15; stype = FE_PENTA15G8; }
		else if (strcmp(sztype, "HEX20G8"     ) == 0) { eshape = ET_HEX20; stype = FE_HEX20G8; }
		else if (strcmp(sztype, "QUAD4G8"     ) == 0) { eshape = ET_QUAD4; stype = FE_SHELL_QUAD4G8; }
		else if (strcmp(sztype, "QUAD4G12"    ) == 0) { eshape = ET_QUAD4; stype = FE_SHELL_QUAD4G12; }
		else if (strcmp(sztype, "QUAD8G18"    ) == 0) { eshape = ET_QUAD8; stype = FE_SHELL_QUAD8G18; }
		else if (strcmp(sztype, "QUAD8G27"    ) == 0) { eshape = ET_QUAD8; stype = FE_SHELL_QUAD8G27; }
		else if (strcmp(sztype, "TRI3G6"      ) == 0) { eshape = ET_TRI3; stype = FE_SHELL_TRI3G6; }
		else if (strcmp(sztype, "TRI3G9"      ) == 0) { eshape = ET_TRI3; stype = FE_SHELL_TRI3G9; }
		else if (strcmp(sztype, "TRI6G14"     ) == 0) { eshape = ET_TRI6; stype = FE_SHELL_TRI6G14; }
		else if (strcmp(sztype, "TRI6G21"     ) == 0) { eshape = ET_TRI6; stype = FE_SHELL_TRI6G21; }
		else
		{
			assert(false);
		}
	}

	// this is a hack to choose between 2D elements and shell elements.
	// NOTE: This is only used by quad/tri elements.
	// TODO: find a better way
	int NDIM = 3;
	if (GetModuleName() == "fluid") NDIM = 2;

	// determine the element type
	FE_Element_Type etype = FE_ELEM_INVALID_TYPE;
	switch (eshape)
	{
	case ET_HEX8   : etype = m_nhex8; break;
	case ET_PENTA6 : etype = FE_PENTA6G6; break;
	case ET_PENTA15: etype = FE_PENTA15G21; break;
	case ET_PYRA5  : etype = FE_PYRA5G8; break;
	case ET_TET4   : etype = m_ntet4; break;
	case ET_TET10  : etype = m_ntet10; break;
	case ET_TET15  : etype = m_ntet15; break;
	case ET_TET20  : etype = m_ntet20; break;
	case ET_HEX20  : etype = (stype == FE_HEX20G8 ? FE_HEX20G8 : FE_HEX20G27); break;
	case ET_HEX27  : etype = FE_HEX27G27; break;
	case ET_QUAD4  : etype = (NDIM == 3 ? stype : FE2D_QUAD4G4); break;
	case ET_TRI3   : etype = (NDIM == 3 ? stype : FE2D_TRI3G1); break;
	case ET_TRI6   : etype = (NDIM == 3 ? stype : FE2D_TRI6G3); break;
	case ET_QUAD8  : etype = (NDIM == 3 ? stype : FE2D_QUAD8G9); break;
	case ET_QUAD9  : etype = FE2D_QUAD9G9; break;
	case ET_TRUSS2 : etype = FE_TRUSS; break;
	default:
		assert(false);
	}

	// determine the element class
	FE_Element_Class eclass = FEElementLibrary::GetElementClass(etype);

	// return the spec
	FE_Element_Spec spec;
	spec.eclass = eclass;
	spec.eshape = eshape;
	spec.etype = etype;
	spec.m_bthree_field_hex = m_b3field_hex;
	spec.m_bthree_field_tet = m_b3field_tet;
    spec.m_bthree_field_shell = m_b3field_shell;
	spec.m_but4 = m_but4;
	spec.m_shell_formulation = m_default_shell;

	// Make sure this is a valid element specification
	assert(FEElementLibrary::IsValid(spec));

	return spec;
}
