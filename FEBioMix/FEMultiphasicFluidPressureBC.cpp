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
#include "FEMultiphasicFluidPressureBC.h"
#include "FEBioMix.h"
#include "FEMultiphasic.h"
#include <FECore/FEModel.h>
#include <FECore/FESurface.h>

//=============================================================================
BEGIN_FECORE_CLASS(FEMultiphasicFluidPressureBC, FEPrescribedSurface)
    ADD_PARAMETER(m_p, "pressure")->setUnits("P")->setLongName("fluid pressure");
    ADD_PARAMETER(m_bshellb , "shell_bottom");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FEMultiphasicFluidPressureBC::FEMultiphasicFluidPressureBC(FEModel* pfem) : FEPrescribedSurface(pfem)
{
    m_p = 0;
    m_dofP = -1;
    m_dofC = -1;
    m_Rgas = 0;
    m_Tabs = 0;
    m_bshellb = false;
}

//-----------------------------------------------------------------------------
//! initialize
bool FEMultiphasicFluidPressureBC::Init()
{
    FEModel* pfem = GetFEModel();
    m_dofP = m_bshellb ? pfem->GetDOFIndex("q") : pfem->GetDOFIndex("p");
    
    m_dofC = m_bshellb ? pfem->GetDOFIndex("shell concentration", 0) : pfem->GetDOFIndex("concentration", 0);
    SetDOFList(m_dofP);

    if (FEPrescribedSurface::Init() == false) return false;
    m_Rgas = GetFEModel()->GetGlobalConstant("R");
    m_Tabs = GetFEModel()->GetGlobalConstant("T");

    m_pe.assign(GetSurface()->Nodes(), 0.0);
    
    return true;
}

//-----------------------------------------------------------------------------
//! Evaluate and prescribe the resistance pressure
void FEMultiphasicFluidPressureBC::Update()
{
    // prescribe this dilatation at the nodes
    FESurface* ps = GetSurface();

    int N = ps->Nodes();
    std::vector<vector<double>> peNodes(N, vector<double>());

    // Project sum of all kappa and osm_coef values from int points to nodes on surface
    // All values put into map, including duplicates
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
        FEElement* e = el.m_elem[0].pe;
        FEMaterial* pm = GetFEModel()->GetMaterial(e->GetMatID());
        FESoluteInterface* psi = pm->ExtractProperty<FESoluteInterface>();
        FESolidElement* se = dynamic_cast<FESolidElement*>(e);
        FEShellElement* sh = dynamic_cast<FEShellElement*>(e);
        if (se) {
            double peo[FEElement::MAX_NODES] = {0};
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
                    // evaluate effective fluid pressure at this node
                    peo[j] = p - m_Rgas*m_Tabs*osc*osm;
                }
            }
            else {
                // loop over face nodes
                for (int j=0; j<el.Nodes(); ++j) {
                    FENode& node = ps->Node(el.m_lnode[j]);
                    // evaluate effective fluid pressure at this node
                    peo[j] = p;
                }
            }
            // only keep the dilatations at the nodes of the surface face
            for (int j=0; j<el.Nodes(); ++j)
                peNodes[el.m_lnode[j]].push_back(peo[j]);
        }
        //If no solid element, insert all 0s
        else if (sh) {
            double peo[FEElement::MAX_NODES] = {0};
            if (psi) {
                const int nsol = psi->Solutes();
                std::vector<double> kappa(nsol,0);
                double osc = 0;
                const int nint = sh->GaussPoints();
                // get the average osmotic coefficient and partition coefficients in the solid element
                for (int j=0; j<nint; ++j) {
                    FEMaterialPoint* pt = sh->GetMaterialPoint(j);
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
                    // evaluate effective fluid pressure at this node
                    peo[j] = p - m_Rgas*m_Tabs*osc*osm;
                }
            }
            else {
                // loop over face nodes
                for (int j=0; j<el.Nodes(); ++j) {
                    FENode& node = ps->Node(el.m_lnode[j]);
                    // evaluate effective fluid pressure at this node
                    peo[j] = p;
                }
            }
            // only keep the dilatations at the nodes of the surface face
            for (int j=0; j<el.Nodes(); ++j)
                peNodes[el.m_lnode[j]].push_back(peo[j]);
        }
    }
    
    //For each node, average the nodal pe
    for (int i=0; i<ps->Nodes(); ++i)
    {
        double pe = 0;
        for (int j = 0; j < peNodes[i].size(); ++j)
            pe += peNodes[i][j];
        pe /= peNodes[i].size();
            
        // store value for now
        m_pe[i] = pe;
    }
 
    FEPrescribedSurface::Update();
}

//-----------------------------------------------------------------------------
void FEMultiphasicFluidPressureBC::GetNodalValues(int nodelid, std::vector<double>& val)
{
    val[0] = m_pe[nodelid];
}

//-----------------------------------------------------------------------------
// copy data from another class
void FEMultiphasicFluidPressureBC::CopyFrom(FEBoundaryCondition* pbc)
{
    // TODO: implement this
    assert(false);
}

//-----------------------------------------------------------------------------
//! serialization
void FEMultiphasicFluidPressureBC::Serialize(DumpStream& ar)
{
    FEPrescribedSurface::Serialize(ar);
    ar & m_pe;
    if (ar.IsShallow()) return;
    ar & m_dofC & m_dofP;
    ar & m_Rgas & m_Tabs;
}
