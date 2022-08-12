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
    ADD_PARAMETER(m_p, "pressure");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FEMultiphasicFSIPressureBC::FEMultiphasicFSIPressureBC(FEModel* pfem) : FEPrescribedSurface(pfem)
{
    m_pfs = nullptr;
    m_p = 0;
    m_dofEF = -1;
    m_dofC = -1;
}

//-----------------------------------------------------------------------------
//! initialize
bool FEMultiphasicFSIPressureBC::Init()
{
    m_dofEF = GetDOFIndex(FEBioMultiphasicFSI::GetVariableName(FEBioMultiphasicFSI::FLUID_DILATATION), 0);
    m_dofC = GetDOFIndex(FEBioMultiphasicFSI::GetVariableName(FEBioMultiphasicFSI::FLUID_CONCENTRATION), 0);
    SetDOFList(m_dofEF);

    if (FEPrescribedSurface::Init() == false) return false;

    // get fluid from first surface element
    // assuming the entire surface bounds the same fluid
    FESurfaceElement& el = GetSurface()->Element(0);
    FEElement* pe = el.m_elem[0];
    if (pe == nullptr) return false;

    // get the material
    FEMaterial* pm = GetFEModel()->GetMaterial(pe->GetMatID());
    m_pfs = dynamic_cast<FEMultiphasicFSI*>(pm);
    if (m_pfs == nullptr) return false;

    m_e.assign(GetSurface()->Nodes(), 0.0);

    return true;
}

//-----------------------------------------------------------------------------
//! Evaluate and prescribe the resistance pressure
void FEMultiphasicFSIPressureBC::Update()
{
    // prescribe this dilatation at the nodes
    FESurface* ps = GetSurface();
    int nsol = 0;
    double T = 0;
    double R = 0;

    nsol = m_pfs->Solutes();
    T = m_pfs->m_Tabs;
    R = m_pfs->m_Rgas;

    int N = ps->Nodes();
    std::vector<vector<double>> oscNodes(N, vector<double>());
    std::vector<vector<double>> caNodes(N, vector<double>());

    //Project sum of all ca and osc values from int points to nodes on surface
    //All values put into map, including duplicates
    for (int i = 0; i < ps->Elements(); ++i)
    {
        FESurfaceElement& el = ps->Element(i);
        FEElement* e = el.m_elem[0];
        FESolidElement* se = dynamic_cast<FESolidElement*>(e);
        if (se) {
            double osci[FEElement::MAX_INTPOINTS];
            double osco[FEElement::MAX_NODES];
            double cai[FEElement::MAX_INTPOINTS];
            double cao[FEElement::MAX_NODES];
            for (int j = 0; j < se->GaussPoints(); ++j) {
                FEMaterialPoint* pt = se->GetMaterialPoint(j);
                FEMultiphasicFSIMaterialPoint* fsp = pt->ExtractData<FEMultiphasicFSIMaterialPoint>();
                if (fsp)
                {
                    osci[j] = m_pfs->GetOsmoticCoefficient()->OsmoticCoefficient(*pt);
                    cai[j] = fsp->m_ca[0];
                    for (int isol = 1; isol < nsol; ++isol)
                        cai[j] += fsp->m_ca[isol];
                }
                else
                {
                    osci[j] = 0;
                    cai[j] = 0;
                }
            }
            // project stresses from integration points to nodes
            se->project_to_nodes(osci, osco);
            se->project_to_nodes(cai, cao);
            // only keep the stresses at the nodes of the contact face
            for (int j = 0; j < el.Nodes(); ++j)
            {
                oscNodes[el.m_lnode[j]].push_back(osco[j]);
                caNodes[el.m_lnode[j]].push_back(cao[j]);
            }
        }
        //If no solid element, insert all 0s
        else {
            for (int j = 0; j < el.Nodes(); ++j)
            {
                oscNodes[el.m_lnode[j]].push_back(0);
                caNodes[el.m_lnode[j]].push_back(0);
            }
        }
    }
    //For each node average the nodal ca and osc and then calculate ef based on desired p
    m_e.assign(N, 0.0);
    FEFluid* pfl = m_pfs->Fluid();
    for (int i = 0; i < ps->Nodes(); ++i)
    {
        //get osmotic component of pressure
        double ca = 0;
        double osc = 0;
        for (int j = 0; j < caNodes[i].size(); ++j)
        {
            ca += caNodes[i][j];
            osc += oscNodes[i][j];
        }
        ca /= caNodes[i].size();
        osc /= caNodes[i].size();

        double c = osc * ca * R;

        //get correct ef for desired pressure
        double e = 0;
        bool good = pfl->Dilatation(T, m_p, c, e);
        assert(good);

        m_e[i] = e;
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
    ar& m_pfs;
    ar& m_dofC;
    ar& m_dofEF;
    ar& m_e;
}
