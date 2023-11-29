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
#include "stdafx.h"
#include "FEMultiphasicFSIPressureBC.h"
#include "FEMultiphasicFSI.h"
#include "FEBioMultiphasicFSI.h"
#include <FECore/FEModel.h>

//=============================================================================
BEGIN_FECORE_CLASS(FEMultiphasicFSIPressureBC, FEPrescribedSurface)
    ADD_PARAMETER(m_p, "pressure")->setUnits("P")->setLongName("fluid pressure");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FEMultiphasicFSIPressureBC::FEMultiphasicFSIPressureBC(FEModel* pfem) : FEPrescribedSurface(pfem)
{
    m_p = 0;
    m_dofEF = -1;
    m_dofC = -1;
    m_Rgas = 0;
    m_Tabs = 0;
}

//-----------------------------------------------------------------------------
//! initialize
bool FEMultiphasicFSIPressureBC::Init()
{
    m_dofEF = GetDOFIndex(FEBioMultiphasicFSI::GetVariableName(FEBioMultiphasicFSI::FLUID_DILATATION), 0);
    m_dofC = GetDOFIndex(FEBioMultiphasicFSI::GetVariableName(FEBioMultiphasicFSI::FLUID_CONCENTRATION), 0);
    SetDOFList(m_dofEF);

    if (FEPrescribedSurface::Init() == false) return false;

    m_Rgas = GetFEModel()->GetGlobalConstant("R");
    m_Tabs = GetFEModel()->GetGlobalConstant("T");
    
    m_e.assign(GetSurface()->Nodes(), 0.0);

    return true;
}

//-----------------------------------------------------------------------------
//! Evaluate and prescribe the resistance pressure
void FEMultiphasicFSIPressureBC::Update()
{
    // prescribe this dilatation at the nodes
    FESurface* ps = GetSurface();
    
    int N = ps->Nodes();
    std::vector<vector<double>> efNodes(N, vector<double>());
    std::vector<vector<double>> caNodes(N, vector<double>());
    
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
        FEElement* e = el.m_elem[0];
        FEMaterial* pm = GetFEModel()->GetMaterial(e->GetMatID());
        FEFluid* pfl = pm->ExtractProperty<FEFluid>();
        FESoluteInterface* psi = pm->ExtractProperty<FESoluteInterface>();
        FESolidElement* se = dynamic_cast<FESolidElement*>(e);
        if (se) {
            double efo[FEElement::MAX_NODES] = {0};
            if (psi) {
                const int nsol = psi->Solutes();
                std::vector<double> kappa(nsol,0);
                double osc = 0;
                const int nint = se->GaussPoints();
                // get the average osmotic coefficient and partition coefficients in the solid element
                for (int j=0; j<nint; ++j) {
                    FEMaterialPoint* pt = se->GetMaterialPoint(j);
                    osc += psi->GetOsmoticCoefficient()->OsmoticCoefficient(*pt);
                    for (int k=0; k<nsol; ++k)
                        kappa[k] += psi->GetPartitionCoefficient(*pt, k);
                }
                osc /= nint;
                for (int k=0; k<nsol; ++k) kappa[k] /= nint;
                // loop over face nodes
                for (int j=0; j<el.Nodes(); ++j) {
                    double osm = 0;
                    FENode& node = ps->Node(el.m_lnode[j]);
                    // calculate osmolarity at this node, using nodal effective solute concentrations
                    for (int k=0; k<nsol; ++k)
                        osm += node.get(m_dofC+psi->GetSolute(k)->GetSoluteID()-1)*kappa[k];
                    // evaluate dilatation at this node
                    bool good = pfl->Dilatation(0, p - m_Rgas*m_Tabs*osc*osm, efo[j]);
                    assert(good);
                }
            }
            else {
                // loop over face nodes
                for (int j=0; j<el.Nodes(); ++j) {
                    FENode& node = ps->Node(el.m_lnode[j]);
                    // evaluate dilatation at this node
                    bool good = pfl->Dilatation(0, p, efo[j]);
                    assert(good);
                }
            }
            // only keep the dilatations at the nodes of the surface face
            for (int j=0; j<el.Nodes(); ++j)
                efNodes[el.m_lnode[j]].push_back(efo[j]);
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
// return the value for node i, dof j
void FEMultiphasicFSIPressureBC::GetNodalValues(int nodelid, std::vector<double>& val)
{
    val[0] = m_e[nodelid];
}

//-----------------------------------------------------------------------------
// copy data from another class
void FEMultiphasicFSIPressureBC::CopyFrom(FEBoundaryCondition* pbc)
{
    // TODO: implement this
    assert(false);
}

//-----------------------------------------------------------------------------
//! serialization
void FEMultiphasicFSIPressureBC::Serialize(DumpStream& ar)
{
    FEPrescribedSurface::Serialize(ar);
    if (ar.IsLoading()) {
        m_e.assign(GetSurface()->Nodes(), 0.0);
    }
    ar & m_e;
    if (ar.IsShallow()) return;
    ar & m_Rgas & m_Tabs;
    ar & m_dofEF & m_dofC;
}
