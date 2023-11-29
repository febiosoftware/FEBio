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
#include "FEInitialFluidPressureTemperature.h"
#include "FEBioThermoFluid.h"
#include "FEThermoFluid.h"
#include "FEFluid.h"
#include <FECore/FEModel.h>
#include <FECore/FEAnalysis.h>

//=============================================================================
BEGIN_FECORE_CLASS(FEInitialFluidPressureTemperature, FEInitialCondition)
// material properties
    ADD_PARAMETER(m_Pdata, "pressure"   )->setUnits(UNIT_PRESSURE);
    ADD_PARAMETER(m_Tdata, "temperature")->setUnits(UNIT_RELATIVE_TEMPERATURE);
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEInitialFluidPressureTemperature::FEInitialFluidPressureTemperature(FEModel* fem) : FENodalIC(fem)
{
}

//-----------------------------------------------------------------------------
bool FEInitialFluidPressureTemperature::Init()
{
    if (SetPDOF("ef") == false) return false;
	if (SetTDOF("T") == false) return false;
    m_e.assign(m_nodeSet->Size(), 0.0);
    m_T.assign(m_nodeSet->Size(), 0.0);

    FEDofList dofs(GetFEModel());
    if (dofs.AddVariable(FEBioThermoFluid::GetVariableName(FEBioThermoFluid::FLUID_DILATATION)) == false) return false;
    if (dofs.AddVariable(FEBioThermoFluid::GetVariableName(FEBioThermoFluid::TEMPERATURE)) == false) return false;
    SetDOFList(dofs);
    
	return FENodalIC::Init();
}

//-----------------------------------------------------------------------------
void FEInitialFluidPressureTemperature::SetTDOF(int ndof) { m_dofT = ndof; }

//-----------------------------------------------------------------------------
bool FEInitialFluidPressureTemperature::SetTDOF(const char* szdof)
{
    FEModel* fem = GetFEModel();
    int ndof = fem->GetDOFIndex(szdof);
    assert(ndof >= 0);
    if (ndof < 0) return false;
    SetTDOF(ndof);
    return true;
}

//-----------------------------------------------------------------------------
void FEInitialFluidPressureTemperature::SetPDOF(int ndof) { m_dofEF = ndof; }

//-----------------------------------------------------------------------------
bool FEInitialFluidPressureTemperature::SetPDOF(const char* szdof)
{
    FEModel* fem = GetFEModel();
    int ndof = fem->GetDOFIndex(szdof);
    assert(ndof >= 0);
    if (ndof < 0) return false;
    SetPDOF(ndof);
    return true;
}

//-----------------------------------------------------------------------------
void FEInitialFluidPressureTemperature::Activate()
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
        // get material
        FEMaterial* pm = GetFEModel()->GetMaterial(el.GetMatID());
        FEThermoFluid* ptfl = pm->ExtractProperty<FEThermoFluid>();
        FEFluid* pfl = pm->ExtractProperty<FEFluid>();
        if (ptfl) {
            double efi[FEElement::MAX_INTPOINTS] = {0};
            double efo[FEElement::MAX_NODES] = {0};
            bool good = true;
            for (int j=0; j<el.GaussPoints(); ++j) {
                FEMaterialPoint* pt = el.GetMaterialPoint(j);
                good = good && ptfl->Dilatation(m_Tdata(*pt), m_Pdata(*pt), efi[j]);
            }
            // project dilatations from integration points to nodes
            el.project_to_nodes(efi, efo);
            if (good) {
                for (int j=0; j<el.Nodes(); ++j)
                    efNodes[el.m_node[j]].push_back(efo[j]);
            }
            else {
                for (int j=0; j<el.Nodes(); ++j) {
                    efo[j] = 0;
                    efNodes[el.m_node[j]].push_back(efo[j]);
                }
            }
        }
        else if (pfl) {
            double efi[FEElement::MAX_INTPOINTS] = {0};
            double efo[FEElement::MAX_NODES] = {0};
            bool good = true;
            for (int j=0; j<el.GaussPoints(); ++j) {
                FEMaterialPoint* pt = el.GetMaterialPoint(j);
                good = good && pfl->Dilatation(0, m_Pdata(*pt), efi[j]);
            }
            // project dilatations from integration points to nodes
            el.project_to_nodes(efi, efo);
            if (good) {
                for (int j=0; j<el.Nodes(); ++j)
                    efNodes[el.m_node[j]].push_back(efo[j]);
            }
            else {
                for (int j=0; j<el.Nodes(); ++j) {
                    efo[j] = 0;
                    efNodes[el.m_node[j]].push_back(efo[j]);
                }
            }
        }
        else break;
    }
    
    //For each node, average the nodal ef
    for (int i=0; i<nodeList.Size(); ++i)
    {
        double ef = 0;
        for (int j = 0; j < efNodes[i].size(); ++j)
            ef += efNodes[i][j];
        ef /= efNodes[i].size();
            
        // store value for now
        m_e[i] = ef;
        
        // get the nodal temperature
        FENode& node = *nodeList.Node(i);
        FEMaterialPoint mp;
        mp.m_r0 = node.m_r0;
        mp.m_index = i;

        m_T[i] = m_Tdata(mp);
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
void FEInitialFluidPressureTemperature::GetNodalValues(int inode, std::vector<double>& values)
{
    values[0] = m_e[inode];
    values[1] = m_T[inode];
}

//-----------------------------------------------------------------------------
void FEInitialFluidPressureTemperature::Serialize(DumpStream& ar)
{
    FENodalIC::Serialize(ar);
    if (ar.IsLoading()) {
        m_e.assign(m_nodeSet->Size(), 0.0);
        m_T.assign(m_nodeSet->Size(), 0.0);
    }
    ar & m_e & m_T;
    if (ar.IsShallow()) return;
    ar & m_dofEF & m_dofT;
}

