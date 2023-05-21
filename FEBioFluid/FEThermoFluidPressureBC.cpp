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
#include "FEThermoFluidPressureBC.h"
#include "FEBioThermoFluid.h"
#include <FECore/FEAnalysis.h>
#include <FECore/FEModel.h>
#include <FECore/FESurface.h>
#include "FEThermoFluid.h"
#include "FEFluid.h"

//=============================================================================
BEGIN_FECORE_CLASS(FEThermoFluidPressureBC, FEPrescribedSurface)
    ADD_PARAMETER(m_p, "pressure")->setUnits("P")->setLongName("fluid pressure");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FEThermoFluidPressureBC::FEThermoFluidPressureBC(FEModel* pfem) : FEPrescribedSurface(pfem)
{
    m_p = 0;
    m_dofEF = -1;
    m_dofT = -1;
    m_Rgas = 0;
    m_Tabs = 0;
    m_psurf = nullptr;
}

//-----------------------------------------------------------------------------
//! initialize
bool FEThermoFluidPressureBC::Init()
{
    m_dofEF = GetDOFIndex(FEBioThermoFluid::GetVariableName(FEBioThermoFluid::FLUID_DILATATION), 0);
    m_dofT = GetDOFIndex(FEBioThermoFluid::GetVariableName(FEBioThermoFluid::TEMPERATURE), 0);
    SetDOFList(m_dofEF);

    if (FEPrescribedSurface::Init() == false) return false;
    m_Rgas = GetFEModel()->GetGlobalConstant("R");
    m_Tabs = GetFEModel()->GetGlobalConstant("T");
    
    m_psurf = FESurfaceBC::GetSurface();

    m_e.assign(GetSurface()->Nodes(), 0.0);

    return true;
}

//-----------------------------------------------------------------------------
//! Evaluate and prescribe the resistance pressure
void FEThermoFluidPressureBC::Update()
{
    // prescribe this dilatation at the nodes
    FESurface* ps = GetSurface();
    FEModel* fem = GetFEModel();
    FEMesh& mesh = fem->GetMesh();
    FETimeInfo& tp = fem->GetTime();

    int N = ps->Nodes();
    std::vector<vector<double>> efNodes(N, vector<double>());

    //Project sum of all ca and osc values from int points to nodes on surface
    //All values put into map, including duplicates
    for (int i=0; i<ps->Elements(); ++i)
    {
        FESurfaceElement& el = ps->Element(i);
        // evaluate average prescribed pressure on this face
        double p = 0;
        for (int j=0; j<el.GaussPoints(); ++j) {
            FEMaterialPoint* pt = el.GetMaterialPoint(j);
            p += m_p(*pt);
        }
        p /= el.GaussPoints();
        // get surface underlying material
        FEElement* e = el.m_elem[0];
        FESolidElement* se = dynamic_cast<FESolidElement*>(e);
        if (se) {
            FEMaterial* pm = GetFEModel()->GetMaterial(e->GetMatID());
            FEThermoFluid* ptfl = pm->ExtractProperty<FEThermoFluid>();
            FEFluid* pfl = pm->ExtractProperty<FEFluid>();
            if (ptfl) {
                // evaluate average temperature on this face
                double T = 0;
                for (int j=0; j<el.Nodes(); ++j) {
                    FENode& node = mesh.Node(el.m_node[j]);
                    T += node.get(m_dofT)*tp.alphaf + node.get_prev(m_dofT)*(1-tp.alphaf);
                }
                T /= el.Nodes();
                double efi[FEElement::MAX_INTPOINTS] = {0};
                double efo[FEElement::MAX_NODES] = {0};
/*                bool good = true;
                for (int j=0; j<se->GaussPoints(); ++j) {
                    FEMaterialPoint* pt = se->GetMaterialPoint(j);
                    good = good && ptfl->Dilatation(T, p, 0, efi[j]);
                    if (!good) break;
                }
                // only keep the dilatations at the nodes of the surface face
                if (good) {
                    // project dilatations from integration points to nodes
                    se->project_to_nodes(efi, efo);
                    for (int j=0; j<el.Nodes(); ++j)
                        efNodes[el.m_lnode[j]].push_back(efo[j]);
                }
                else {*/
                    for (int j=0; j<el.Nodes(); ++j) {
                        FENode& node = mesh.Node(el.m_node[j]);
                        efo[j] = node.get_prev(m_dofEF);
                        ptfl->Dilatation(T, p, efo[j]);
                        efNodes[el.m_lnode[j]].push_back(efo[j]);
                    }
//                }
            }
            else if (pfl) {
                double efi[FEElement::MAX_INTPOINTS] = {0};
                double efo[FEElement::MAX_NODES] = {0};
                bool good = true;
                for (int j=0; j<se->GaussPoints(); ++j) {
                    FEMaterialPoint* pt = se->GetMaterialPoint(j);
                    good = good && pfl->Dilatation(0, p, efi[j]);
                    if (!good) break;
                }
                // only keep the dilatations at the nodes of the surface face
                if (good) {
                    // project dilatations from integration points to nodes
                    se->project_to_nodes(efi, efo);
                    for (int j=0; j<el.Nodes(); ++j)
                        efNodes[el.m_lnode[j]].push_back(efo[j]);
                }
                else {
                    for (int j=0; j<el.Nodes(); ++j) {
                        FENode& node = mesh.Node(el.m_node[j]);
                        efo[j] = node.get_prev(m_dofEF);
                        pfl->Dilatation(0, p, efo[j]);
                        efNodes[el.m_lnode[j]].push_back(efo[j]);
                    }
                }
            }
            else break;
        }
        //If no solid element, insert all 0s
        else {
            for (int j=0; j<el.Nodes(); ++j)
                efNodes[el.m_lnode[j]].push_back(0);
        }
    }
    
    //For each node, average the nodal ef
    for (int i=0; i<ps->Nodes(); ++i)
    {
        double ef = 0;
        for (int j = 0; j < efNodes[i].size(); ++j)
            ef += efNodes[i][j];
        ef /= efNodes[i].size();
            
        // store value for now
        m_e[i] = ef;
    }
 
    FEPrescribedSurface::Update();
}

//-----------------------------------------------------------------------------
void FEThermoFluidPressureBC::GetNodalValues(int nodelid, std::vector<double>& val)
{
    val[0] = m_e[nodelid];
}

//-----------------------------------------------------------------------------
// copy data from another class
void FEThermoFluidPressureBC::CopyFrom(FEBoundaryCondition* pbc)
{
    // TODO: implement this
    assert(false);
}

//-----------------------------------------------------------------------------
//! serialization
void FEThermoFluidPressureBC::Serialize(DumpStream& ar)
{
    FEPrescribedSurface::Serialize(ar);
    ar & m_e;
    if (ar.IsShallow()) return;
    ar & m_dofT & m_dofEF;
    ar & m_Rgas & m_Tabs;
    ar & m_psurf;
}
