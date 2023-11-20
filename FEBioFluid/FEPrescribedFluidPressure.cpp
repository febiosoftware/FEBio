/*This file is part of the FEBio source code and is licensed under the MIT license
 listed below.

 See Copyright-FEBio.txt for details.

 Copyright (c) 2020 University of Utah, The Trustees of Columbia University in
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
#include "FEPrescribedFluidPressure.h"
#include "FEBioFluid.h"
#include <FECore/FESurface.h>
#include <FECore/FEModel.h>
#include "FEFluid.h"

 //-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEPrescribedFluidPressure, FEPrescribedSurface)
    ADD_PARAMETER(m_p, "pressure")->setUnits(UNIT_PRESSURE)->setLongName("fluid pressure");
END_FECORE_CLASS();


//-----------------------------------------------------------------------------
 //! constructor
FEPrescribedFluidPressure::FEPrescribedFluidPressure(FEModel* fem) : FEPrescribedSurface(fem)
{
	m_p = 0.0;
    m_psurf = nullptr;
}

//-----------------------------------------------------------------------------
//! initialize
bool FEPrescribedFluidPressure::Init()
{
    // get the dof index
    m_dofEF = GetDOFIndex(FEBioFluid::GetVariableName(FEBioFluid::FLUID_DILATATION), 0);
    if (m_dofEF < 0) return false;
    SetDOFList(m_dofEF);

    if (FEPrescribedSurface::Init() == false) return false;
    
    m_psurf = GetSurface();

    m_e.assign(GetSurface()->Nodes(), 0.0);
    
    // do an initial Update so that the dilatations are set properly at the very first time step
    Update();
    
    return true;
}

//-----------------------------------------------------------------------------
void FEPrescribedFluidPressure::UpdateDilatation()
{
	int N = m_psurf->Nodes();
	std::vector<vector<double>> efNodes(N, vector<double>());

	//Project sum of all ca and osc values from int points to nodes on surface
	//All values put into map, including duplicates
	for (int i = 0; i < m_psurf->Elements(); ++i)
	{
		FESurfaceElement& el = m_psurf->Element(i);
		// evaluate average prescribed pressure on this face
		double p = 0;
		for (int j = 0; j < el.GaussPoints(); ++j) {
			FEMaterialPoint* pt = el.GetMaterialPoint(j);
			p += m_p(*pt);
		}
		p /= el.GaussPoints();
		FEElement* e = el.m_elem[0];
		FEMaterial* pm = GetFEModel()->GetMaterial(e->GetMatID());
		FEFluid* pfl = pm->ExtractProperty<FEFluid>();
		if (pfl == nullptr) break;
		FESolidElement* se = dynamic_cast<FESolidElement*>(e);
		if (se) {
			double efi[FEElement::MAX_INTPOINTS] = { 0 };
			double efo[FEElement::MAX_NODES] = { 0 };
			for (int j = 0; j < se->GaussPoints(); ++j) {
				FEMaterialPoint* pt = se->GetMaterialPoint(j);
				bool good = pfl->Dilatation(0, p, efi[j]);
				assert(good);
			}
			// project dilatations from integration points to nodes
			se->project_to_nodes(efi, efo);
			// only keep the dilatations at the nodes of the surface face
			for (int j = 0; j < el.Nodes(); ++j)
				efNodes[el.m_lnode[j]].push_back(efo[j]);
		}
		//If no solid element, insert all 0s
		else {
			for (int j = 0; j < el.Nodes(); ++j)
				efNodes[el.m_lnode[j]].push_back(0);
		}
	}
	//For each node, average the nodal ef
	for (int i = 0; i < m_psurf->Nodes(); ++i)
	{
		double ef = 0;
		for (int j = 0; j < efNodes[i].size(); ++j)
			ef += efNodes[i][j];
		ef /= efNodes[i].size();

		// store value for now
		m_e[i] = ef;
	}
}

//-----------------------------------------------------------------------------
//! Evaluate and prescribe the resistance pressure
void FEPrescribedFluidPressure::UpdateModel() { Update(); }
void FEPrescribedFluidPressure::Update()
{
	UpdateDilatation();
    FEPrescribedSurface::Update();

	// TODO: Is this necessary?
	GetFEModel()->SetMeshUpdateFlag(true);
}

//-----------------------------------------------------------------------------
void FEPrescribedFluidPressure::PrepStep(std::vector<double>& ui, bool brel)
{
	UpdateDilatation();
	FEPrescribedSurface::PrepStep(ui, brel);
}

//-----------------------------------------------------------------------------
void FEPrescribedFluidPressure::GetNodalValues(int nodelid, std::vector<double>& val)
{
    val[0] = m_e[nodelid];

	FENode& node = GetMesh().Node(m_nodeList[nodelid]);
	node.set(m_dofEF, m_e[nodelid]);

}

void FEPrescribedFluidPressure::CopyFrom(FEBoundaryCondition* pbc)
{
    // TODO: implement
    assert(false);
}

//-----------------------------------------------------------------------------
//! serialization
void FEPrescribedFluidPressure::Serialize(DumpStream& ar)
{
    FEPrescribedSurface::Serialize(ar);
    ar & m_e;
    if (ar.IsShallow()) return;
    ar & m_dofEF;
    ar & m_psurf;
}
