/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2023 University of Utah, The Trustees of Columbia University in
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
#include "FEFluidSolutesPressureLC.h"
#include "FEBioFluidSolutes.h"
#include "FEFluidMaterial.h"
#include "FESolutesMaterial.h"
#include <FECore/FECoreKernel.h>
#include <FECore/FEModel.h>
#include <FECore/FEMesh.h>
#include <FECore/FEMaterial.h>
#include <FEBioMix/FESoluteInterface.h>

BEGIN_FECORE_CLASS(FEFluidSolutesPressureLC, FESurfaceConstraint)
    ADD_PARAMETER(m_lc.m_rhs, "pressure");
    ADD_PARAMETER(m_lc.m_laugon, "laugon");
    ADD_PARAMETER(m_lc.m_tol, "tol");
    ADD_PARAMETER(m_lc.m_eps, "penalty");
    ADD_PARAMETER(m_lc.m_naugmin, "minaug");
    ADD_PARAMETER(m_lc.m_naugmax, "maxaug");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEFluidSolutesPressureLC::FEFluidSolutesPressureLC(FEModel* pfem) : FESurfaceConstraint(pfem), m_surf(pfem), m_lc(pfem)
{
}

//-----------------------------------------------------------------------------
//! Initializes data structures.
void FEFluidSolutesPressureLC::Activate()
{
    // don't forget to call base class
    FENLConstraint::Activate();
    m_lc.Activate();
}

//-----------------------------------------------------------------------------
bool FEFluidSolutesPressureLC::Init()
{
    // initialize surface
    m_surf.Init();
    FEModel* fem = GetFEModel();
    m_dofEF = fem->GetDOFIndex(FEBioFluidSolutes::GetVariableName(FEBioFluidSolutes::FLUID_DILATATION), 0);
    m_dofC = fem->GetDOFIndex(FEBioFluidSolutes::GetVariableName(FEBioFluidSolutes::FLUID_CONCENTRATION), 0);
    
    // assume all faces have same underlying material as first face of surface
    // get surface element
    FESurfaceElement& el = m_surf.Element(0);
    // get underlying solid element
    FESolidElement* pe = dynamic_cast<FESolidElement*>(el.m_elem[0]);
    if (pe == nullptr) return false;
    FEMaterial* pm = GetFEModel()->GetMaterial(pe->GetMatID());
    FEFluidMaterial* pf = pm->ExtractProperty<FEFluidMaterial>();
    // get the local solute id
    FESoluteInterface* psi = dynamic_cast<FESoluteInterface*>(pm);
    if (psi == nullptr) return false;
    // get fluid bulk modulus
    double K = pf->BulkModulus(*pe->GetMaterialPoint(0));
    double RTPhi = fem->GetGlobalConstant("R")*fem->GetGlobalConstant("T");
    int nsol = psi->Solutes();

    int N = m_surf.Nodes();
    
    // create linear constraints
    for (int i=0; i<N; ++i) {
        
        //  -K ef - R T Phi kappa ctilde = p
        FEAugLagLinearConstraint* pLC0 = fecore_alloc(FEAugLagLinearConstraint, fem);
        FENode& node = m_surf.Node(i);
        pLC0->AddDOF(node.GetID(), m_dofEF, -1.0);
        for (int j=0; j<nsol; ++j) {
            int sid = psi->GetSolute(j)->GetSoluteID() - 1;
            pLC0->AddDOF(node.GetID(), m_dofC + sid, RTPhi/K);
        }
        // add the linear constraint to the system
        m_lc.add(pLC0);
    }
    m_lc.m_rhs /= K;
    
    return true;
}

//-----------------------------------------------------------------------------
void FEFluidSolutesPressureLC::Serialize(DumpStream& ar)
{
    m_lc.Serialize(ar);
    if (ar.IsShallow()) return;
    ar & m_dofEF & m_dofC;
}

//-----------------------------------------------------------------------------
void FEFluidSolutesPressureLC::LoadVector(FEGlobalVector& R, const FETimeInfo& tp) { m_lc.LoadVector(R, tp); }
void FEFluidSolutesPressureLC::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp) { m_lc.StiffnessMatrix(LS, tp); }
bool FEFluidSolutesPressureLC::Augment(int naug, const FETimeInfo& tp) { return m_lc.Augment(naug, tp); }
void FEFluidSolutesPressureLC::BuildMatrixProfile(FEGlobalMatrix& M) { m_lc.BuildMatrixProfile(M); }
