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
#include "FEThermoFluidPressureLoad.h"
#include "FEBioThermoFluid.h"
#include <FECore/FEModel.h>
#include <FECore/FELinearSystem.h>

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEThermoFluidPressureLoad, FESurfaceConstraint)
    ADD_PARAMETER(m_p, "pressure")->setUnits("P")->setLongName("fluid pressure");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEThermoFluidPressureLoad::FEThermoFluidPressureLoad(FEModel* pfem) : FESurfaceConstraint(pfem), m_surf(pfem)
{
    m_pfluid = nullptr;
    m_p = 0;
    
    m_dofEF = pfem->GetDOFIndex(FEBioThermoFluid::GetVariableName(FEBioThermoFluid::FLUID_DILATATION), 0);
    m_dofT = pfem->GetDOFIndex(FEBioThermoFluid::GetVariableName(FEBioThermoFluid::TEMPERATURE), 0);
}

//-----------------------------------------------------------------------------
bool FEThermoFluidPressureLoad::Init()
{
    if (FESurfaceConstraint::Init() == false) return false;
    
    m_surf.Init();
    
    // get fluid from first surface element
    // assuming the entire surface bounds the same fluid
    FESurfaceElement& el = m_surf.Element(0);
    FEElement* pe = el.m_elem[0];
    if (pe == nullptr) return false;
    
    // get the material
    FEMaterial* pm = GetFEModel()->GetMaterial(pe->GetMatID());
    m_pfluid = pm->ExtractProperty<FEThermoFluid>();
    if (m_pfluid == nullptr) return false;
    
    m_EQ.resize(m_surf.Nodes(), -1);
    m_Lm.resize(m_surf.Nodes(), 0.0);
    m_Lmp.resize(m_surf.Nodes(), 0.0);

    return true;
}

//-----------------------------------------------------------------------------
// allocate equations
int FEThermoFluidPressureLoad::InitEquations(int neq)
{
    int n = neq;
    for (int i = 0; i < m_surf.Nodes(); ++i) m_EQ[i] = n++;

    return n - neq;
}

//-----------------------------------------------------------------------------
void FEThermoFluidPressureLoad::UnpackLM(vector<int>& lm, int n)
{
    int ndof = 3;
    
    // add the dofs of rigid body A
    lm.reserve(ndof);
    FENode& node = m_surf.Node(n);
    lm.push_back(node.m_ID[m_dofEF]);
    lm.push_back(node.m_ID[m_dofT]);
    lm.push_back(m_EQ[n]);
}

//-----------------------------------------------------------------------------
// Build the matrix profile
void FEThermoFluidPressureLoad::BuildMatrixProfile(FEGlobalMatrix& M)
{
    for (int i=0; i<m_surf.Nodes(); ++i) {
        vector<int> lm;
        UnpackLM(lm, i);
        
        // add it to the pile
        M.build_add(lm);
    }
}

//-----------------------------------------------------------------------------
void FEThermoFluidPressureLoad::Update(const std::vector<double>& Ui, const std::vector<double>& ui)
{
    for (int i = 0; i < m_surf.Nodes(); ++i)
    {
        if (m_EQ[i] != -1) m_Lm[i] = m_Lmp[i] + Ui[m_EQ[i]] + ui[m_EQ[i]];
    }
}

void FEThermoFluidPressureLoad::PrepStep()
{
    for (int i = 0; i < m_surf.Nodes(); ++i)
    {
        m_Lmp[i] = m_Lm[i];
    }
}

void FEThermoFluidPressureLoad::UpdateIncrements(std::vector<double>& Ui, const std::vector<double>& ui)
{
    for (int i = 0; i < m_surf.Nodes(); ++i)
    {
        if (m_EQ[i] != -1) Ui[m_EQ[i]] += ui[m_EQ[i]];
    }
}

//-----------------------------------------------------------------------------
//! Evaluate and prescribe the resistance pressure
void FEThermoFluidPressureLoad::Update()
{
}

//-----------------------------------------------------------------------------
//! serialization
void FEThermoFluidPressureLoad::Serialize(DumpStream& ar)
{
    FESurfaceConstraint::Serialize(ar);
    if (ar.IsShallow()) return;
    ar & m_pfluid;
    ar & m_Lm & m_Lmp;
    ar & m_dofT & m_dofEF;
}

//-----------------------------------------------------------------------------
//! \todo Why is this class not using the FESolver for assembly?
void FEThermoFluidPressureLoad::LoadVector(FEGlobalVector& R, const FETimeInfo& tp)
{
    int ndof = 3;
    vector<double> fe(ndof, 0.0);
    
    // for now assume that the prescribed pressure is the same on the entire surface
    FESurfaceElement& el = m_surf.Element(0);
    // evaluate average prescribed pressure on this face
    double p0 = 0;
    for (int j=0; j<el.GaussPoints(); ++j) {
        FEMaterialPoint* pt = el.GetMaterialPoint(j);
        p0 += m_p(*pt);
    }
    p0 /= el.GaussPoints();

    double alpha = tp.alphaf;
    
    for (int i=0; i<m_surf.Nodes(); ++i) {
        FENode& node = m_surf.Node(i);
        double e = node.get(m_dofEF)*alpha + node.get_prev(m_dofEF)*(1-alpha);
        double T = node.get(m_dofT)*alpha + node.get_prev(m_dofT)*(1-alpha);
        double p = m_pfluid->GetElastic()->Pressure(e, T);
        double dpJ = m_pfluid->GetElastic()->Tangent_Strain(e, T);
        double dpT = m_pfluid->GetElastic()->Tangent_Temperature(e, T);
        double lam = m_Lm[i]*alpha + m_Lmp[i]*(1-alpha);
        fe[0] = lam*dpJ;
        fe[1] = lam*dpT;
        fe[2] = p - p0;
        vector<int> lm;
        UnpackLM(lm,i);
        R.Assemble(lm, fe);
    }
    
}

//-----------------------------------------------------------------------------
//! \todo Why is this class not using the FESolver for assembly?
void FEThermoFluidPressureLoad::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
{
    int ndof = 3;
    FEElementMatrix ke;
    ke.resize(ndof, ndof);
    double alpha = tp.alphaf;
    
    for (int i=0; i<m_surf.Nodes(); ++i) {
        ke.zero();
        FENode& node = m_surf.Node(i);
        double e = node.get(m_dofEF)*alpha + node.get_prev(m_dofEF)*(1-alpha);
        double T = node.get(m_dofT)*alpha + node.get_prev(m_dofT)*(1-alpha);
        double dpJ = m_pfluid->GetElastic()->Tangent_Strain(e, T);
        double dpT = m_pfluid->GetElastic()->Tangent_Temperature(e, T);
        double lam = m_Lm[i]*alpha + m_Lmp[i]*(1-alpha);
        double dpJ2 = m_pfluid->GetElastic()->Tangent_Strain_Strain(e, T);
        double dpJT = m_pfluid->GetElastic()->Tangent_Strain_Temperature(e, T);
        double dpT2 = m_pfluid->GetElastic()->Tangent_Temperature_Temperature(e, T);
        
        mat3d Kab(lam*dpJ2, lam*dpJT, dpJ,
                  lam*dpJT, lam*dpT2, dpT,
                  dpJ, dpT, 0);
        ke.sub(0, 0, Kab);

        // unpack LM
        vector<int> lm;
        UnpackLM(lm, i);
        ke.SetIndices(lm);
        
        // assemle into global stiffness matrix
        LS.Assemble(ke);
    }
}
