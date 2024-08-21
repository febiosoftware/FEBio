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
#include "FEThermoFluidSolver.h"
#include <FECore/log.h>
#include <FECore/FEModel.h>
#include <FECore/FELinearSystem.h>

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEThermoFluidPressureLoad, FESurfaceConstraint)
    ADD_PARAMETER(m_p, "pressure")->setUnits("P")->setLongName("fluid pressure");
    ADD_PARAMETER(m_laugon, "laugon"        )->setLongName("Enforcement method")->setEnums("PENALTY\0AUGLAG\0LAGMULT\0");
    BEGIN_PARAM_GROUP("Augmentation");
        ADD_PARAMETER(m_tol, "tol")->setLongName("tolerance");
        ADD_PARAMETER(m_eps, "penalty");
        ADD_PARAMETER(m_naugmin, "minaug");
        ADD_PARAMETER(m_naugmax, "maxaug");
    END_PARAM_GROUP();
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEThermoFluidPressureLoad::FEThermoFluidPressureLoad(FEModel* pfem) : FESurfaceConstraint(pfem), m_surf(pfem)
{
    m_pfluid = nullptr;
    m_p = 0;
    m_tol = 0.1;
    m_laugon = FECore::PENALTY_METHOD;
    m_eps = 1.0;
    m_naugmin = 0;
    m_naugmax = 10;
    
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
    FEElement* pe = el.m_elem[0].pe;
    if (pe == nullptr) return false;
    
    // get the material
    FEMaterial* pm = GetFEModel()->GetMaterial(pe->GetMatID());
    m_pfluid = pm->ExtractProperty<FEThermoFluid>();
    if (m_pfluid == nullptr) return false;
    
    switch (m_laugon) {
        case 0:
        case 1:
            m_Lm.resize(m_surf.Nodes(), 0.0);
            m_Lmp.resize(m_surf.Nodes(), 0.0);
            break;
        case 2:
        {
            m_EQ.resize(m_surf.Nodes(), -1);
            m_Lm.resize(m_surf.Nodes(), 0.0);
            m_Lmp.resize(m_surf.Nodes(), 0.0);
        }
            break;
        default:
            break;
    }
    
    // project the prescribed pressure to the nodes
    m_surf.ProjectToNodes(m_p, m_pn);
    
    return true;
}

//-----------------------------------------------------------------------------
// allocate equations
int FEThermoFluidPressureLoad::InitEquations(int neq)
{
    if (m_laugon < FECore::LAGMULT_METHOD) return 0;
    int n = neq;
    for (int i = 0; i < m_surf.Nodes(); ++i) m_EQ[i] = n++;

    return n - neq;
}

//-----------------------------------------------------------------------------
void FEThermoFluidPressureLoad::UnpackLM(vector<int>& lm, int n)
{
    int ndof;
    if (m_laugon < FECore::LAGMULT_METHOD) ndof = 2;
    else ndof = 3;
    
    lm.reserve(ndof);
    FENode& node = m_surf.Node(n);
    lm.push_back(node.m_ID[m_dofEF]);
    lm.push_back(node.m_ID[m_dofT]);
    if (m_laugon == FECore::LAGMULT_METHOD) lm.push_back(m_EQ[n]);
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
    if (m_laugon < FECore::LAGMULT_METHOD) return;
    for (int i = 0; i < m_surf.Nodes(); ++i)
    {
        if (m_EQ[i] != -1) m_Lm[i] = m_Lmp[i] + Ui[m_EQ[i]] + ui[m_EQ[i]];
    }
}

//-----------------------------------------------------------------------------
void FEThermoFluidPressureLoad::PrepStep()
{
    m_surf.ProjectToNodes(m_p, m_pn);
    if (m_laugon < FECore::LAGMULT_METHOD) return;
    for (int i = 0; i < m_surf.Nodes(); ++i)
    {
        m_Lmp[i] = m_Lm[i];
    }
}

//-----------------------------------------------------------------------------
void FEThermoFluidPressureLoad::UpdateIncrements(std::vector<double>& Ui, const std::vector<double>& ui)
{
    if (m_laugon < FECore::LAGMULT_METHOD) return;
    for (int i = 0; i < m_surf.Nodes(); ++i)
    {
        if (m_EQ[i] != -1) Ui[m_EQ[i]] += ui[m_EQ[i]];
    }
}

//-----------------------------------------------------------------------------
//! serialization
void FEThermoFluidPressureLoad::Serialize(DumpStream& ar)
{
    FESurfaceConstraint::Serialize(ar);
    if (ar.IsShallow()) return;
    ar & m_pfluid;
    ar & m_dofT & m_dofEF;
    if (m_laugon == FECore::AUGLAG_METHOD) ar & m_Lm;
    else if (m_laugon > FECore::AUGLAG_METHOD)
        ar & m_Lm & m_Lmp;
}

//-----------------------------------------------------------------------------
//! \todo Why is this class not using the FESolver for assembly?
void FEThermoFluidPressureLoad::LoadVector(FEGlobalVector& R, const FETimeInfo& tp)
{
    double alpha = tp.alphaf;
    
    if (m_laugon == FECore::LAGMULT_METHOD) {
        int ndof = 3;
        vector<double> fe(ndof, 0.0);
        
        for (int i=0; i<m_surf.Nodes(); ++i) {
            FENode& node = m_surf.Node(i);
            double e = node.get(m_dofEF)*alpha + node.get_prev(m_dofEF)*(1-alpha);
            double T = node.get(m_dofT)*alpha + node.get_prev(m_dofT)*(1-alpha);
            double p = m_pfluid->GetElastic()->Pressure(e, T);
            double dpJ = m_pfluid->GetElastic()->Tangent_Strain(e, T);
            double dpT = m_pfluid->GetElastic()->Tangent_Temperature(e, T);
            double lam = m_Lm[i]*alpha + m_Lmp[i]*(1-alpha);
            double f = p - m_pn[i];
            fe[0] = -lam*dpJ;
            fe[1] = -lam*dpT;
            fe[2] = -f;
            vector<int> lm;
            UnpackLM(lm,i);
            R.Assemble(lm, fe);
        }
    }
    else {
        int ndof = 2;
        vector<double> fe(ndof, 0.0);
        
        for (int i=0; i<m_surf.Nodes(); ++i) {
            FENode& node = m_surf.Node(i);
            double e = node.get(m_dofEF)*alpha + node.get_prev(m_dofEF)*(1-alpha);
            double T = node.get(m_dofT)*alpha + node.get_prev(m_dofT)*(1-alpha);
            double p = m_pfluid->GetElastic()->Pressure(e, T);
            double dpJ = m_pfluid->GetElastic()->Tangent_Strain(e, T);
            double dpT = m_pfluid->GetElastic()->Tangent_Temperature(e, T);
            double p0 = m_pn[i];
            double c = m_Lm[i] + m_eps*(p - p0);
            fe[0] = -c*dpJ;
            fe[1] = -c*dpT;
            vector<int> lm;
            UnpackLM(lm,i);
            R.Assemble(lm, fe);
        }
    }

}

//-----------------------------------------------------------------------------
//! \todo Why is this class not using the FESolver for assembly?
void FEThermoFluidPressureLoad::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
{
    double alpha = tp.alphaf;
    
    if (m_laugon == FECore::LAGMULT_METHOD) {
        int ndof = 3;
        FEElementMatrix ke;
        ke.resize(ndof, ndof);
        for (int i=0; i<m_surf.Nodes(); ++i) {
            ke.zero();
            FENode& node = m_surf.Node(i);
            double e = node.get(m_dofEF)*alpha + node.get_prev(m_dofEF)*(1-alpha);
            double T = node.get(m_dofT)*alpha + node.get_prev(m_dofT)*(1-alpha);
            double p = m_pfluid->GetElastic()->Pressure(e, T);
            double dpJ = m_pfluid->GetElastic()->Tangent_Strain(e, T);
            double dpT = m_pfluid->GetElastic()->Tangent_Temperature(e, T);
            double lam = m_Lm[i]*alpha + m_Lmp[i]*(1-alpha);
            double f = p - m_pn[i];
            double dpJ2 = m_pfluid->GetElastic()->Tangent_Strain_Strain(e, T);
            double dpJT = m_pfluid->GetElastic()->Tangent_Strain_Temperature(e, T);
            double dpT2 = m_pfluid->GetElastic()->Tangent_Temperature_Temperature(e, T);
            
            mat3d Kab(lam*dpJ2, lam*dpJT, dpJ,
                      lam*dpJT, lam*dpT2, dpT,
                      dpJ, dpT, 0);
            ke.add(0, 0, Kab);
            
            // unpack LM
            vector<int> lm;
            UnpackLM(lm, i);
            ke.SetIndices(lm);
            
            // assemle into global stiffness matrix
            LS.Assemble(ke);
        }
    }
    else {
        int ndof = 2;
        FEElementMatrix ke;
        ke.resize(ndof, ndof);
        for (int i=0; i<m_surf.Nodes(); ++i) {
            ke.zero();
            FENode& node = m_surf.Node(i);
            double e = node.get(m_dofEF)*alpha + node.get_prev(m_dofEF)*(1-alpha);
            double T = node.get(m_dofT)*alpha + node.get_prev(m_dofT)*(1-alpha);
            double p = m_pfluid->GetElastic()->Pressure(e, T);
            double p0 = m_pn[i];
            double dpJ = m_pfluid->GetElastic()->Tangent_Strain(e, T);
            double dpT = m_pfluid->GetElastic()->Tangent_Temperature(e, T);
            double c = m_Lm[i] + m_eps*(p - p0);
            double dpJ2 = m_pfluid->GetElastic()->Tangent_Strain_Strain(e, T);
            double dpJT = m_pfluid->GetElastic()->Tangent_Strain_Temperature(e, T);
            double dpT2 = m_pfluid->GetElastic()->Tangent_Temperature_Temperature(e, T);
            
            ke(0, 0) += m_eps*dpJ*dpJ + c*dpJ2;
            ke(0, 1) += m_eps*dpJ*dpT + c*dpJT;
            ke(1, 0) += m_eps*dpJ*dpT + c*dpJT;
            ke(1, 1) += m_eps*dpT*dpT + c*dpT2;

            // unpack LM
            vector<int> lm;
            UnpackLM(lm, i);
            ke.SetIndices(lm);
            
            // assemle into global stiffness matrix
            LS.Assemble(ke);
        }
    }
}

//-----------------------------------------------------------------------------
bool FEThermoFluidPressureLoad::Augment(int naug, const FETimeInfo& tp)
{
    if (m_laugon != FECore::AUGLAG_METHOD) return true;
    
    double alpha = tp.alphaf;
    
    // calculate lag multipliers
    double L0 = 0, L1 = 0;
    for (int i=0; i<m_surf.Nodes(); ++i) {
        FENode& node = m_surf.Node(i);
        double e = node.get(m_dofEF)*alpha + node.get_prev(m_dofEF)*(1-alpha);
        double T = node.get(m_dofT)*alpha + node.get_prev(m_dofT)*(1-alpha);
        double p = m_pfluid->GetElastic()->Pressure(e, T);
        double lam = m_Lm[i] + m_eps*(p - m_pn[i]);
        L0 += m_Lm[i]*m_Lm[i];
        L1 += lam*lam;
        m_Lmp[i] = lam;
    }
    
    L0 = sqrt(L0);
    L1 = sqrt(L1);
    
    double d;
    if (L1 != 0)
        d = fabs((L1 - L0)/L1);
    else d = fabs(L1 - L0);
    
    const std::string name = GetName();
    feLog("constraint %s: %15.7lg %15.7lg %15.7lg\n", name.c_str(), L0, fabs(L1 - L0), fabs(m_tol*L1));
    
    bool bconv = false;
    if (d <= m_tol) bconv = true;
    if ((m_naugmax >= 0) && (naug >= m_naugmax)) bconv = true;
    if (naug < m_naugmin) bconv = false;
    
    if (bconv == false)
    {
        for (int i=0; i<m_surf.Nodes(); ++i) m_Lm[i] = m_Lmp[i];
    }
    
    return bconv;
}

