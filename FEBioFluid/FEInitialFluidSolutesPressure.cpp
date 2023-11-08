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
#include "FEInitialFluidSolutesPressure.h"
#include "FEBioFluidSolutes.h"
#include "FEFluidSolutes.h"
#include <FECore/FEModel.h>

//=============================================================================
BEGIN_FECORE_CLASS(FEInitialFluidSolutesPressure, FEInitialCondition)
// material properties
    ADD_PARAMETER(m_Pdata, "value"   )->setUnits(UNIT_PRESSURE);
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEInitialFluidSolutesPressure::FEInitialFluidSolutesPressure(FEModel* fem) : FENodalIC(fem)
{
    m_dofEF = -1;
    m_dofC = -1;
    m_Rgas = 0;
    m_Tabs = 0;
    m_Fc = 0;
}

//-----------------------------------------------------------------------------
bool FEInitialFluidSolutesPressure::Init()
{
    if (SetPDOF("ef") == false) return false;
    m_e.assign(m_nodeSet->Size(), 0.0);

    FEDofList dofs(GetFEModel());
    if (dofs.AddVariable(FEBioFluidSolutes::GetVariableName(FEBioFluidSolutes::FLUID_DILATATION)) == false) return false;
    SetDOFList(dofs);
    m_dofC = GetDOFIndex(FEBioFluidSolutes::GetVariableName(FEBioFluidSolutes::FLUID_CONCENTRATION), 0);
    m_Rgas = GetFEModel()->GetGlobalConstant("R");
    m_Tabs = GetFEModel()->GetGlobalConstant("T");
    m_Fc   = GetFEModel()->GetGlobalConstant("Fc");

	return FENodalIC::Init();
}

//-----------------------------------------------------------------------------
void FEInitialFluidSolutesPressure::SetPDOF(int ndof) { m_dofEF = ndof; }

//-----------------------------------------------------------------------------
bool FEInitialFluidSolutesPressure::SetPDOF(const char* szdof)
{
    FEModel* fem = GetFEModel();
    int ndof = fem->GetDOFIndex(szdof);
    assert(ndof >= 0);
    if (ndof < 0) return false;
    SetPDOF(ndof);
    return true;
}

//-----------------------------------------------------------------------------
void FEInitialFluidSolutesPressure::Activate()
{
    // prescribe this dilatation at the nodes
    FEModel* fem = GetFEModel();
    
    FENodeList nodeList = m_nodeSet->GetNodeList();

    int N = nodeList.Size();
    std::vector<vector<double>> efNodes(N, vector<double>());
    
    FEMesh& mesh = *m_nodeSet->GetMesh();

    // evaluate average prescribed pressure and temperature in each element
    // project them from int points to nodes
    for (int i=0; i<mesh.Elements(); ++i)
    {
        FEElement& el = *mesh.Element(i);
        FESolidElement* se = dynamic_cast<FESolidElement*>(&el);
        // get material
        FEMaterial* pm = GetFEModel()->GetMaterial(el.GetMatID());
        FEFluid* pfl = pm->ExtractProperty<FEFluid>();
        FESoluteInterface* psi = pm->ExtractProperty<FESoluteInterface>();
        if (pfl) {
            double efi[FEElement::MAX_INTPOINTS] = {0};
            double efo[FEElement::MAX_NODES] = {0};
            if (psi) {
                const int nsol = psi->Solutes();
                std::vector<double> kappa(nsol,0);
                double osc = 0, p = 0, epot = 0;
                const int nint = se->GaussPoints();
                // get the average osmotic coefficient and partition coefficients in the solid element
                for (int j=0; j<nint; ++j) {
                    FEMaterialPoint* pt = se->GetMaterialPoint(j);
                    osc += psi->GetOsmoticCoefficient()->OsmoticCoefficient(*pt);
                    p += m_Pdata(*pt);
                    epot += psi->GetElectricPotential(*pt);
                    for (int k=0; k<nsol; ++k)
                        kappa[k] += psi->GetSolute(k)->m_pSolub->Solubility(*pt);
                }
                osc /= nint;
                p /= nint;
                epot /= nint;
                double ex = exp(-m_Fc*epot/m_Rgas/m_Tabs);
                for (int k=0; k<nsol; ++k) kappa[k] *= pow(ex,psi->GetSolute(k)->ChargeNumber())/nint;
                // loop over element nodes
                for (int j=0; j<el.Nodes(); ++j) {
                    double osm = 0;
                    FENode& node = mesh.Node(el.m_lnode[j]);
                    // calculate osmolarity at this node, using nodal effective solute concentrations
                    for (int k=0; k<nsol; ++k)
                        osm += node.get(m_dofC+psi->GetSolute(k)->GetSoluteID()-1)*kappa[k];
                    bool good = pfl->Dilatation(0, p - m_Rgas*m_Tabs*osc*osm, efo[j]);
                    assert(good);
                }
            }
            else {
                double p = 0;
                const int nint = se->GaussPoints();
                // get the average osmotic coefficient and partition coefficients in the solid element
                for (int j=0; j<nint; ++j) {
                    FEMaterialPoint* pt = se->GetMaterialPoint(j);
                    p += m_Pdata(*pt);
                }
                p /= nint;
                // loop over element nodes
                for (int j=0; j<el.Nodes(); ++j) {
                    FENode& node = mesh.Node(el.m_lnode[j]);
                    bool good = pfl->Dilatation(0, p, efo[j]);
                    assert(good);
                }
            }
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
    for (int i=0; i<mesh.Nodes(); ++i)
    {
        double ef = 0;
        for (int j = 0; j < efNodes[i].size(); ++j)
            ef += efNodes[i][j];
        ef /= efNodes[i].size();
        
        // store value for now
        m_e[i] = ef;
    }

    FEStepComponent::Activate();
    if (m_dofs.IsEmpty()) return;

    int dofs = (int)m_dofs.Size();
    std::vector<double> val(dofs, 0.0);

    for (int i = 0; i<N; ++i)
    {
        FENode& node = *m_nodeSet->Node(i);

        // get the nodal values
        GetNodalValues(i, val);
        
        for (int j = 0; j < dofs; ++j)
        {
            node.set(m_dofs[j], val[j]);
        }
    }
}

//-----------------------------------------------------------------------------
void FEInitialFluidSolutesPressure::GetNodalValues(int inode, std::vector<double>& values)
{
    values[0] = m_e[inode];
}

//-----------------------------------------------------------------------------
void FEInitialFluidSolutesPressure::Serialize(DumpStream& ar)
{
    FENodalIC::Serialize(ar);
    if (ar.IsLoading())
        m_e.assign(m_nodeSet->Size(), 0.0);
    ar & m_e;
    if (ar.IsShallow()) return;
    ar & m_dofEF;
    ar & m_dofC;
    ar & m_Rgas & m_Tabs & m_Fc;
}

