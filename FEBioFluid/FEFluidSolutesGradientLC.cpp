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
#include "FEFluidSolutesGradientLC.h"
#include "FEBioFluidSolutes.h"
#include "FEFluidMaterial.h"
#include "FESolutesMaterial.h"
#include "FEFluidSolutesDomain3D.h"
#include <FECore/FECoreKernel.h>
#include <FECore/FEModel.h>
#include <FECore/FEMesh.h>
#include <FECore/FEMaterial.h>
#include <FEBioMix/FESoluteInterface.h>

BEGIN_FECORE_CLASS(FEFluidSolutesGradientLC, FESurfaceConstraint)
    ADD_PARAMETER(m_isol   , "solute_id")->setEnums("$(solutes)");
    ADD_PARAMETER(m_lc.m_laugon, "laugon");
    ADD_PARAMETER(m_lc.m_tol, "tol");
    ADD_PARAMETER(m_lc.m_eps, "penalty");
    ADD_PARAMETER(m_lc.m_naugmin, "minaug");
    ADD_PARAMETER(m_lc.m_naugmax, "maxaug");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEFluidSolutesGradientLC::FEFluidSolutesGradientLC(FEModel* pfem) : FESurfaceConstraint(pfem), m_surf(pfem), m_lc(pfem)
{
    m_isol = -1;
    m_binit = false;
}

//-----------------------------------------------------------------------------
//! Initializes data structures.
void FEFluidSolutesGradientLC::Activate()
{
    // don't forget to call base class
    FESurfaceConstraint::Activate();
    
    if (m_binit == false)
    {
        FEModel* fem = GetFEModel();
        m_dofC = fem->GetDOFIndex(FEBioFluidSolutes::GetVariableName(FEBioFluidSolutes::FLUID_CONCENTRATION), 0);
        
        int NE = m_surf.Elements();
        
        // create linear constraints
        for (int i=0; i<NE; ++i) {
            // get surface element
            FESurfaceElement& el = m_surf.Element(i);
            // get surface normal
            vec3d nu(0,0,0);
            for (int n=0; n<el.GaussPoints(); ++n) {
                FESurfaceMaterialPoint* pt = dynamic_cast<FESurfaceMaterialPoint*>(el.GetMaterialPoint(n));
                nu += pt->dxr ^ pt->dxs;
            }
            nu.unit();
            
            // get underlying solid element
            FESolidElement* pe = dynamic_cast<FESolidElement*>(el.m_elem[0]);
            if (pe == nullptr) break;
            // determine the solid domain to which this solid element belongs
            FESolidDomain* sdom = dynamic_cast<FESolidDomain*>(pe->GetMeshPartition());
            // get element data
            int nint = pe->GaussPoints();
            int neln = pe->Nodes();
            
            FEMaterial* pm = GetFEModel()->GetMaterial(pe->GetMatID());
            // make sure that we have solutes in this selement
            FESoluteInterface* psi = dynamic_cast<FESoluteInterface*>(pm);
            if (psi == nullptr) break;
            
            // for each integration point in the solid element
            for (int n=0; n<nint; ++n) {
                FEMaterialPoint& pt = *pe->GetMaterialPoint(n);
                FEFluidMaterialPoint& pb = *(pt.ExtractData<FEFluidMaterialPoint>());
                FEFluidSolutesMaterialPoint& ps = *(pt.ExtractData<FEFluidSolutesMaterialPoint>());
                
                // get shape function derivatives
                double* Gr = pe->Gr(n);
                double* Gs = pe->Gs(n);
                double* Gt = pe->Gt(n);
                
                // get contravariant basis vectors
                vec3d gcntv[3];
                sdom->ContraBaseVectors(*pe, n, gcntv);
                
                // evaluate gradient of shape function
                vector<vec3d> gradN(neln);
                for (int i=0; i<neln; ++i)
                    gradN[i] = gcntv[0]*Gr[i] + gcntv[1]*Gs[i] + gcntv[2]*Gt[i];
                
                //  sum c gradN.nu = 0
                FEAugLagLinearConstraint* pLC0 = fecore_alloc(FEAugLagLinearConstraint, fem);
                int sid = m_isol - 1;
                for (int i=0; i<neln; ++i) {
                    FENode node = sdom->Node(pe->m_lnode[i]);
                    pLC0->AddDOF(node.GetID(), m_dofC + sid, gradN[i]*nu);
                }
                // add the linear constraint to the system
                m_lc.add(pLC0);
            }
        }
        m_lc.m_rhs = 0;
        
        m_lc.Init();
        m_lc.Activate();
        
        m_binit = true;
    }
}

//-----------------------------------------------------------------------------
bool FEFluidSolutesGradientLC::Init()
{
    if (m_isol <= 0) return false;
    
    // initialize surface
    return m_surf.Init();
}

//-----------------------------------------------------------------------------
void FEFluidSolutesGradientLC::Serialize(DumpStream& ar)
{
    m_lc.Serialize(ar);
    if (ar.IsShallow()) return;
    ar & m_dofC;
}

//-----------------------------------------------------------------------------
void FEFluidSolutesGradientLC::LoadVector(FEGlobalVector& R, const FETimeInfo& tp) { m_lc.LoadVector(R, tp); }
void FEFluidSolutesGradientLC::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp) { m_lc.StiffnessMatrix(LS, tp); }
bool FEFluidSolutesGradientLC::Augment(int naug, const FETimeInfo& tp) { return m_lc.Augment(naug, tp); }
void FEFluidSolutesGradientLC::BuildMatrixProfile(FEGlobalMatrix& M) { m_lc.BuildMatrixProfile(M); }
