/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

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
#include <FECore/FELoadCurve.h>
#include <FECore/FESurfaceMap.h>
#include <FECore/FEDomainMap.h>
#include <FECore/FEEdge.h>
#include <FECore/FEConstValueVec3.h>
#include <FECore/log.h>
#include <FECore/FEDataGenerator.h>
#include <FECore/FEModule.h>
#include <FECore/FEPointFunction.h>
#include <FECore/FEBodyLoad.h>
#include <FECore/FESurfacePairConstraint.h>
#include <FECore/FENLConstraint.h>
#include <sstream>

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

	m_nquad4 = FE_QUAD4G4;
	m_nquad8 = FE_QUAD8G9;
	m_nquad9 = FE_QUAD9G9;

	// 3-field formulation flags
	m_b3field_hex = true;
	m_b3field_tet = false;
    m_b3field_shell = false;
    m_b3field_quad = true;
    m_b3field_tri = false;

	// shell formulation
	m_default_shell = NEW_SHELL;
    m_shell_norm_nodal = true;

	// UT4 formulation off by default
	m_but4 = false;
	m_ut4_alpha = 0.05;
	m_ut4_bdev = false;

	// UDG hourglass parameter
	m_udghex_hg = 1.0;
}

//-----------------------------------------------------------------------------
FEModelBuilder::~FEModelBuilder() 
{

}

//-----------------------------------------------------------------------------
void FEModelBuilder::SetActiveModule(const std::string& moduleName)
{
	m_fem.SetActiveModule(moduleName);
}

//-----------------------------------------------------------------------------
//! Get the module name
std::string FEModelBuilder::GetModuleName() const
{
	return m_fem.GetModuleName();
}

//-----------------------------------------------------------------------------
FEAnalysis* FEModelBuilder::CreateNewStep(bool allocSolver)
{
	// default analysis type should match module name
	std::string modName = GetModuleName();
	FEAnalysis* pstep = fecore_new<FEAnalysis>(modName.c_str(), &m_fem);

	// make sure we have a solver defined
	FESolver* psolver = pstep->GetFESolver();
	if ((psolver == 0) && allocSolver)
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
	if (m_defaultSolver.empty() == false) sztype = m_defaultSolver.c_str();
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
//! Get the mesh
FEMesh& FEModelBuilder::GetMesh()
{
	return m_fem.GetMesh();
}

//-----------------------------------------------------------------------------
//! get the FE model
FEModel& FEModelBuilder::GetFEModel()
{
	return m_fem;
}

//-----------------------------------------------------------------------------
//! Create a domain
FEDomain* FEModelBuilder::CreateDomain(FE_Element_Spec espec, FEMaterial* mat)
{
	FECoreKernel& febio = FECoreKernel::GetInstance();
	FEDomain* pdom = febio.CreateDomain(espec, &m_fem.GetMesh(), mat);
	return pdom;
}

//-----------------------------------------------------------------------------
FEAnalysis* FEModelBuilder::GetStep(bool allocSolver)
{
	if (m_pStep == 0)
	{
		m_pStep = CreateNewStep(allocSolver);
		m_fem.AddStep(m_pStep);
		if (m_fem.Steps() == 1)
		{
			m_fem.SetCurrentStep(m_pStep);
			m_fem.SetCurrentStepIndex(0);
		}
	}
	return m_pStep;
}

void FEModelBuilder::AddMaterial(FEMaterial* pmat)
{
	m_fem.AddMaterial(pmat);
}

//-----------------------------------------------------------------------------
void FEModelBuilder::AddComponent(FEStepComponent* pmc)
{
	if (m_nsteps > 0)
	{
		GetStep()->AddStepComponent(pmc);
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
	m_fem.AddModelLoad(pfc);
	AddComponent(pfc);
}

//-----------------------------------------------------------------------------
void FEModelBuilder::AddSurfaceLoad(FESurfaceLoad* psl)
{
	m_fem.AddModelLoad(psl);
	AddComponent(psl);
}

//-----------------------------------------------------------------------------
void FEModelBuilder::AddEdgeLoad(FEEdgeLoad* pel)
{
	m_fem.AddModelLoad(pel);
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
void FEModelBuilder::AddRigidComponent(FEStepComponent* pmc) { assert(false); }

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
	s.Create(fs);

	// read faces
	for (int i = 0; i<faces; ++i)
	{
		FESurfaceElement& el = s.Element(i);
		FEFacetSet::FACET& fi = fs.Face(i);

		// set the element type/integration rule
		if (bnodal)
		{
			if      (fi.ntype == 4) el.SetType(FE_QUAD4NI);
			else if (fi.ntype == 3) el.SetType(FE_TRI3NI);
			else if (fi.ntype == 6) el.SetType(FE_TRI6NI);
			else if (fi.ntype == 8) el.SetType(FE_QUAD8NI);
			else if (fi.ntype == 9) el.SetType(FE_QUAD9NI);
			else return false;
		}
		else
		{
			if      (fi.ntype ==  3) el.SetType(m_ntri3);
			else if (fi.ntype ==  4) el.SetType(m_nquad4);
			else if (fi.ntype ==  6) el.SetType(m_ntri6);
			else if (fi.ntype ==  7) el.SetType(m_ntri7);
			else if (fi.ntype == 10) el.SetType(m_ntri10);
			else if (fi.ntype ==  8) el.SetType(m_nquad8);
			else if (fi.ntype ==  9) el.SetType(m_nquad9);
			else return false;
		}

		int N = el.Nodes(); assert(N == fi.ntype);
		for (int j = 0; j<N; ++j) el.m_node[j] = fi.node[j];
	}

    // copy the name
    s.SetName(fs.GetName());

	// finish building the surface
	s.InitSurface();

	// allocate material point data
	s.CreateMaterialPointData();
    
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
	int varQ = dofs.AddVariable("shell rotation", VAR_VEC3);
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
	int varV = dofs.AddVariable("velocity", VAR_VEC3);
	dofs.SetDOFName(varV, 0, "vx");
	dofs.SetDOFName(varV, 1, "vy");
	dofs.SetDOFName(varV, 2, "vz");
	int varW = dofs.AddVariable("relative fluid velocity", VAR_VEC3);
	dofs.SetDOFName(varW, 0, "wx");
	dofs.SetDOFName(varW, 1, "wy");
	dofs.SetDOFName(varW, 2, "wz");
	int varAW = dofs.AddVariable("relative fluid acceleration", VAR_VEC3);
	dofs.SetDOFName(varAW, 0, "awx");
	dofs.SetDOFName(varAW, 1, "awy");
	dofs.SetDOFName(varAW, 2, "awz");
	int varVF = dofs.AddVariable("fluid velocity", VAR_VEC3);
	dofs.SetDOFName(varVF, 0, "vfx");
	dofs.SetDOFName(varVF, 1, "vfy");
	dofs.SetDOFName(varVF, 2, "vfz");
	int varAF = dofs.AddVariable("fluid acceleration", VAR_VEC3);
	dofs.SetDOFName(varAF, 0, "afx");
	dofs.SetDOFName(varAF, 1, "afy");
	dofs.SetDOFName(varAF, 2, "afz");
	int varEF = dofs.AddVariable("fluid dilatation");
	dofs.SetDOFName(varEF, 0, "ef");
	int varAEF = dofs.AddVariable("fluid dilatation tderiv");
	dofs.SetDOFName(varAEF, 0, "aef");
	int varQV = dofs.AddVariable("shell velocity", VAR_VEC3);
	dofs.SetDOFName(varQV, 0, "svx");
	dofs.SetDOFName(varQV, 1, "svy");
	dofs.SetDOFName(varQV, 2, "svz");
	int varQA = dofs.AddVariable("shell acceleration", VAR_VEC3);
	dofs.SetDOFName(varQA, 0, "sax");
	dofs.SetDOFName(varQA, 1, "say");
	dofs.SetDOFName(varQA, 2, "saz");
    int varT = dofs.AddVariable("temperature");
    dofs.SetDOFName(varT, 0, "T");
	// must be last variable definition!!
	int varC = dofs.AddVariable("concentration", VAR_ARRAY); // we start with zero concentrations
															 // must be last variable definition!!
	int varSC = dofs.AddVariable("shell concentration", VAR_ARRAY); // we start with zero concentrations
    int varAC = dofs.AddVariable("concentration tderiv", VAR_ARRAY); // we start with zero concentrations
    // must be last variable definition!!
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
    else if (strcmp(sztype, "pyra13" ) == 0) eshape = ET_PYRA13;
	else if (strcmp(sztype, "tet4"   ) == 0) eshape = ET_TET4;
	else if (strcmp(sztype, "tet5"   ) == 0) eshape = ET_TET5;
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
	else if (strcmp(sztype, "line2"  ) == 0) eshape = ET_TRUSS2;
	else if (strcmp(sztype, "ut4"    ) == 0) { eshape = ET_TET4; m_but4 = true; }
	else
	{
		// new way for defining element type and integration rule at the same time
		// this is useful for multi-step analyses where the geometry is read in before the control section.
		if      (strcmp(sztype, "TET4G4"      ) == 0) { eshape = ET_TET4 ; m_ntet4  = FE_TET4G4; }
		else if (strcmp(sztype, "TET10G4"     ) == 0) { eshape = ET_TET10; m_ntet10 = FE_TET10G4; }
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
		else if (strcmp(sztype, "HEX8G1"      ) == 0) { eshape = ET_HEX8; m_nhex8 = FE_HEX8G1; }
		else if (strcmp(sztype, "HEX8G8"      ) == 0) { eshape = ET_HEX8; m_nhex8 = FE_HEX8G8; }
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
    case ET_PYRA13 : etype = FE_PYRA13G8; break;
	case ET_TET4   : etype = m_ntet4; break;
	case ET_TET5   : etype = FE_TET5G4; break;
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

	// set three-field flag
	switch (eshape)
	{
	case ET_HEX8   : spec.m_bthree_field = m_b3field_hex; break;
	case ET_PENTA6 : spec.m_bthree_field = m_b3field_hex; break;
	case ET_PYRA5  : spec.m_bthree_field = m_b3field_hex; break;
	case ET_TET4   : spec.m_bthree_field = m_b3field_tet; break;
	case ET_TET5   : spec.m_bthree_field = m_b3field_tet; break;
	case ET_TET10  : spec.m_bthree_field = m_b3field_tet; break;
	case ET_TET15  : spec.m_bthree_field = m_b3field_tet; break;
	case ET_TET20  : spec.m_bthree_field = m_b3field_tet; break;
	case ET_QUAD4  : spec.m_bthree_field = m_b3field_shell || m_b3field_quad; break;
	case ET_TRI3   : spec.m_bthree_field = m_b3field_shell || m_b3field_tri; break;
	case ET_TRI6   : spec.m_bthree_field = m_b3field_shell; break;
	case ET_QUAD8  : spec.m_bthree_field = m_b3field_shell; break;
	case ET_QUAD9  : spec.m_bthree_field = m_b3field_shell; break;
	}

	spec.m_but4 = m_but4;
	spec.m_ut4_alpha = m_ut4_alpha;
	spec.m_ut4_bdev = m_ut4_bdev;
	spec.m_shell_formulation = m_default_shell;
    spec.m_shell_norm_nodal = m_shell_norm_nodal;

	// Make sure this is a valid element specification
	assert(FEElementLibrary::IsValid(spec));

	return spec;
}

void FEModelBuilder::AddMappedParameter(FEParam* p, FECoreBase* parent, const char* szmap, int index)
{
	MappedParameter mp;
	mp.pp = p;
	mp.pc = parent;
	mp.szname = strdup(szmap);
	mp.index = index;

	m_mappedParams.push_back(mp);
}

void FEModelBuilder::AddMeshDataGenerator(FEMeshDataGenerator* gen, FEDomainMap* map, FEParamDouble* pp)
{
	m_mapgen.push_back(DataGen{ gen, map, pp });
}

void FEModelBuilder::ApplyParameterMaps()
{
	FEMesh& mesh = m_fem.GetMesh();
	for (int i = 0; i < m_mappedParams.size(); ++i)
	{
		MappedParameter& mp = m_mappedParams[i];
		FEParam& p = *mp.pp;

		FEDataMap* data = (FEDataMap*)mesh.FindDataMap(mp.szname);
		if (data == nullptr)
		{
			stringstream ss;
			ss << "Can't find map \"" << mp.szname << "\" for parameter \"" << p.name() << "\"";
			throw std::runtime_error(ss.str());
		}

		FEItemList* itemList = data->GetItemList();

		// find the map of this parameter
		if (p.type() == FE_PARAM_DOUBLE_MAPPED)
		{
			FEParamDouble& v = p.value<FEParamDouble>(mp.index);
			FEMappedValue* map = fecore_alloc(FEMappedValue, &m_fem);
			if (data->DataType() != FE_DOUBLE)
			{
				std::stringstream ss;
				ss << "Cannot assign map \"" << data->GetName() << "\" to parameter \"" << p.name() << "\" : bad data type";
				string err = ss.str();
				throw std::runtime_error(err.c_str());
			}
			map->setDataMap(data);
			v.setValuator(map);
			v.SetItemList(itemList);
		}
		else if (p.type() == FE_PARAM_VEC3D_MAPPED)
		{
			FEParamVec3& v = p.value<FEParamVec3>();
			FEMappedValueVec3* map = fecore_alloc(FEMappedValueVec3, &m_fem);
			if (data->DataType() != FE_VEC3D)
			{
				std::stringstream ss;
				ss << "Cannot assign map \"" << data->GetName() << "\" to parameter \"" << p.name() << "\" : bad data type";
				string err = ss.str();
				throw std::runtime_error(err.c_str());
			}
			map->setDataMap(data);
			v.setValuator(map);
			v.SetItemList(itemList);
		}
		else if (p.type() == FE_PARAM_MAT3D_MAPPED)
		{
			FEParamMat3d& v = p.value<FEParamMat3d>();
			FEMappedValueMat3d* map = fecore_alloc(FEMappedValueMat3d, &m_fem);
			if (data->DataType() != FE_MAT3D)
			{
				std::stringstream ss;
				ss << "Cannot assign map \"" << data->GetName() << "\" to parameter \"" << p.name() << "\" : bad data type";
				string err = ss.str();
				throw std::runtime_error(err.c_str());
			}
			map->setDataMap(data);
			v.setValuator(map);
			v.SetItemList(itemList);
		}
        else if (p.type() == FE_PARAM_MAT3DS_MAPPED)
        {
            FEParamMat3ds& v = p.value<FEParamMat3ds>();
            FEMappedValueMat3ds* map = fecore_alloc(FEMappedValueMat3ds, &m_fem);
            if (data->DataType() != FE_MAT3DS)
            {
                std::stringstream ss;
                ss << "Cannot assign map \"" << data->GetName() << "\" to parameter \"" << p.name() << "\" : bad data type";
                string err = ss.str();
                throw std::runtime_error(err.c_str());
            }
            map->setDataMap(data);
            v.setValuator(map);
            v.SetItemList(itemList);
        }
		else
		{
			assert(false);
		}
	}
}

FENodeSet* FEModelBuilder::FindNodeSet(const string& setName)
{
	FEMesh& mesh = m_fem.GetMesh();

	if (setName.compare(0, 9, "@surface:") == 0)
	{
		// see if we can find a surface
		string surfName = setName.substr(9);
		FEFacetSet* surf = mesh.FindFacetSet(surfName);
		if (surf == nullptr) return nullptr;

		// we might have been here before. If so, we already create a nodeset
		// with the same name as the surface, so look for that first.
		FENodeSet* ps = mesh.FindNodeSet(surfName);
		if (ps) return ps;

		// okay, first time here, so let's create a node set from this surface
		FENodeList nodeList = surf->GetNodeList();
		ps = new FENodeSet(&m_fem);
		ps->Add(nodeList);
		ps->SetName(surfName);
		mesh.AddNodeSet(ps);

		return ps;
	}
	if (setName.compare(0, 6, "@edge:") == 0)
	{
		// see if we can find an edge
		string edgeName = setName.substr(6);
		FESegmentSet* edge = mesh.FindSegmentSet(edgeName);
		if (edge == nullptr) return nullptr;

		// we might have been here before. If so, we already create a nodeset
		// with the same name as the edge, so look for that first.
		FENodeSet* ps = mesh.FindNodeSet(edgeName);
		if (ps) return ps;

		// okay, first time here, so let's create a node set from this surface
		FENodeList nodeList = edge->GetNodeList();
		ps = new FENodeSet(&m_fem);
		ps->Add(nodeList);
		ps->SetName(edgeName);
		mesh.AddNodeSet(ps);

		return ps;
	}
	else if (setName.compare(0, 10, "@elem_set:") == 0)
	{
		// see if we can find an element set
		string esetName = setName.substr(10);
		FEElementSet* part = mesh.FindElementSet(esetName);
		if (part == nullptr) return nullptr;

		// we might have been here before. If so, we already create a nodeset
		// with the same name as the surface, so look for that first.
		FENodeSet* ps = mesh.FindNodeSet(esetName);
		if (ps) return ps;

		// okay, first time here, so let's create a node set from this element set
		FENodeList nodeList = part->GetNodeList();
		ps = new FENodeSet(&m_fem);
		ps->Add(nodeList);
		ps->SetName(esetName);
		mesh.AddNodeSet(ps);

		return ps;
	}
	else return mesh.FindNodeSet(setName);
}

void FEModelBuilder::MapLoadCurveToFunction(FEPointFunction* pf, int lc, double scale)
{
	MapLCToFunction m = { lc, scale, pf };
	m_lc2fnc.push_back(m);
}

void FEModelBuilder::ApplyLoadcurvesToFunctions()
{
	FEModel& fem = m_fem;
	for (int i = 0; i < m_lc2fnc.size(); ++i)
	{
		MapLCToFunction& m = m_lc2fnc[i];

		FELoadController* plc = fem.GetLoadController(m.lc); assert(plc);
		FELoadCurve* lc = dynamic_cast<FELoadCurve*>(plc);

		m.pf->SetInterpolation(lc->GetInterpolation());
		m.pf->SetExtendMode(lc->GetExtendMode());
		m.pf->SetPoints(lc->GetPoints());
		if (m.scale != 1.0) m.pf->Scale(m.scale);
	}
}

bool FEModelBuilder::GenerateMeshDataMaps()
{
	FEModel& fem = GetFEModel();
	FEMesh& mesh = GetMesh();
	for (int i = 0; i < m_mapgen.size(); ++i)
	{
		FEMeshDataGenerator* gen = m_mapgen[i].gen;

		// try to initialize the generator
		if (gen->Init() == false) return false;

		FENodeDataGenerator* ngen = dynamic_cast<FENodeDataGenerator*>(gen);
		if (ngen)
		{
			// generate the node data map
			FENodeDataMap* map = ngen->Generate();
			if (map == nullptr) return false;

			map->SetName(ngen->GetName());

			// see if this map is already defined
			string mapName = map->GetName();
			FENodeDataMap* oldMap = dynamic_cast<FENodeDataMap*>(mesh.FindDataMap(mapName));
			if (oldMap)
			{
				// TODO: implement merge
				assert(false);
				// it is, so merge it
//				oldMap->Merge(*map);

				// we can now delete this map
				delete map;
			}
			else
			{
				// nope, so add it
				map->SetName(mapName);
				mesh.AddDataMap(map);
			}
		}

		FEFaceDataGenerator* fgen = dynamic_cast<FEFaceDataGenerator*>(gen);
		if (fgen)
		{
			FESurfaceMap* map = fgen->Generate();
			if (map == nullptr) return false;
			map->SetName(fgen->GetName());
			mesh.AddDataMap(map);
		}

		FEElemDataGenerator* egen = dynamic_cast<FEElemDataGenerator*>(gen);
		if (egen)
		{
			FEDomainMap* map = m_mapgen[i].map;
			FEParamDouble* pp = m_mapgen[i].pp;

			// generate the data
			if (map)
			{
				if (egen->Generate(*map) == false) return false;
			}
			else
			{
				map = egen->Generate();
				map->SetName(egen->GetName());
			}

			// see if this map is already defined
			string mapName = map->GetName();
			FEDomainMap* oldMap = dynamic_cast<FEDomainMap*>(mesh.FindDataMap(mapName));
			if (oldMap)
			{
				// it is, so merge it
				oldMap->Merge(*map);

				// we can now delete this map
				delete map;
			}
			else
			{
				// nope, so add it
				map->SetName(mapName);
				mesh.AddDataMap(map);

				// apply the map
				if (pp)
				{
					FEMappedValue* val = fecore_alloc(FEMappedValue, &fem);
					val->setDataMap(map);
					pp->setValuator(val);
				}
			}
		}
	}

	return true;
}

// finish the build process
bool FEModelBuilder::Finish()
{
	ApplyLoadcurvesToFunctions();
	if (GenerateMeshDataMaps() == false) return false;
	ApplyParameterMaps();
	return true;
}

FEBModel& FEModelBuilder::GetFEBModel()
{
	return m_feb;
}
