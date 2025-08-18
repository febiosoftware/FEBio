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
#include "FETiedFluidInterface.h"
#include <FECore/FENormalProjection.h>
#include <FECore/FEClosestPointProjection.h>
#include <FECore/log.h>
#include <FECore/DumpStream.h>
#include <FECore/FEGlobalMatrix.h>
#include <FECore/FELinearSystem.h>
#include <FECore/FEModel.h>
#include "FEBioFluid.h"

//-----------------------------------------------------------------------------
// FETiedFluidSurface
//-----------------------------------------------------------------------------

FETiedFluidSurface::FETiedFluidSurface(FEModel* pfem) : FESurfaceConstraint(pfem), m_dofWE(pfem), m_surf(pfem), m_lc(pfem)
{
    m_binit = false;
}

//-----------------------------------------------------------------------------
void FETiedFluidSurface::Serialize(DumpStream& ar)
{
//    ar & m_pme;
    ar & m_rs;
}

//-----------------------------------------------------------------------------
void FETiedFluidSurface::Activate()
{
    if (m_binit == false)
    {
        // don't forget to call base class
        FESurfaceConstraint::Activate();
        
        FEModel& fem = *GetFEModel();
        FEMesh& mesh = fem.GetMesh();
        FETiedFluidSurface* ts = GetSibling();
        FESurface* surf = ts->GetSurface();
        FETiedFluidInterface* tfi = GetContactInterface();
        
        // copy settings
        m_lc.m_laugon = tfi->m_laugon;
        m_lc.m_tol = tfi->m_tol;
        m_lc.m_eps = tfi->m_eps;
        m_lc.m_naugmin = tfi->m_naugmin;
        m_lc.m_naugmax = tfi->m_naugmax;

        // create linear constraints
        for (int i=0; i<m_surf.Nodes(); ++i) {
            FENode& node = m_surf.Node(i);
            if (m_pme[i] == nullptr) break;
            if ((node.HasFlags(FENode::EXCLUDE) == false) && (node.m_rid == -1)) {
                
                // get shape functions on sibling surface at projection point
                double Na[FEElement::MAX_NODES];
                m_pme[i]->shape_fnc(Na, m_rs[i][0], m_rs[i][1]);

                // extract the fluid bulk modulus for the solid element underneath this face
                FEElement& el = *(m_pme[i]->m_elem[0].pe);
                FEMaterial* pm = fem.GetMaterial(el.GetMatID());
                FEFluidMaterial* fm = dynamic_cast<FEFluidMaterial*>(pm);
                double K = fm->BulkModulus(*el.GetMaterialPoint(0));

                //  constrain components of velocities on primary and secondary surfaces
                FEAugLagLinearConstraint *pLCvx = fecore_alloc(FEAugLagLinearConstraint, GetFEModel());
                FEAugLagLinearConstraint *pLCvy = fecore_alloc(FEAugLagLinearConstraint, GetFEModel());
                FEAugLagLinearConstraint *pLCvz = fecore_alloc(FEAugLagLinearConstraint, GetFEModel());
                FEAugLagLinearConstraint *pLCef = fecore_alloc(FEAugLagLinearConstraint, GetFEModel());
                pLCvx->AddDOF(node.GetID(), m_dofWE[0], 1);
                pLCvy->AddDOF(node.GetID(), m_dofWE[1], 1);
                pLCvz->AddDOF(node.GetID(), m_dofWE[2], 1);
                pLCef->AddDOF(node.GetID(), m_dofWE[3], 1*K);
                for (int j=0; j<m_pme[i]->Nodes(); ++j) {
                    FENode& n2 = surf->Node(m_pme[i]->m_lnode[j]);
                    pLCvx->AddDOF(n2.GetID(), m_dofWE[0], -Na[j]);
                    pLCvy->AddDOF(n2.GetID(), m_dofWE[1], -Na[j]);
                    pLCvz->AddDOF(n2.GetID(), m_dofWE[2], -Na[j]);
                    pLCef->AddDOF(n2.GetID(), m_dofWE[3], -Na[j]*K);
                }
                m_lc.add(pLCvx);
                m_lc.add(pLCvy);
                m_lc.add(pLCvz);
                m_lc.add(pLCef);
            }
        }
        m_lc.Init();
        m_lc.Activate();
        m_binit = true;
    }
}

//-----------------------------------------------------------------------------
bool FETiedFluidSurface::Init()
{
    // initialize surface data first
    if (FESurfaceConstraint::Init() == false) return false;
    
    FESurface* surf = GetSurface();
    m_pme.assign(surf->Nodes(), nullptr);
    double rs[2] = {0,0};
    m_rs.assign(surf->Nodes(), std::vector<double>(2,0));

	// set the dof list
	if (m_dofWE.AddVariable(FEBioFluid::GetVariableName(FEBioFluid::RELATIVE_FLUID_VELOCITY)) == false) return false;
    if (m_dofWE.AddVariable(FEBioFluid::GetVariableName(FEBioFluid::FLUID_DILATATION)) == false) return false;

    return true;
}

//-----------------------------------------------------------------------------
void FETiedFluidSurface::LoadVector(FEGlobalVector& R, const FETimeInfo& tp) { m_lc.LoadVector(R, tp); }
void FETiedFluidSurface::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp) { m_lc.StiffnessMatrix(LS, tp); }
bool FETiedFluidSurface::Augment(int naug, const FETimeInfo& tp) { return m_lc.Augment(naug, tp); }
void FETiedFluidSurface::BuildMatrixProfile(FEGlobalMatrix& M) { m_lc.BuildMatrixProfile(M); }
int FETiedFluidSurface::InitEquations(int neq) { return m_lc.InitEquations(neq); }
void FETiedFluidSurface::Update(const std::vector<double>& Ui, const std::vector<double>& ui) { m_lc.Update(Ui, ui); }
void FETiedFluidSurface::UpdateIncrements(std::vector<double>& Ui, const std::vector<double>& ui) { m_lc.UpdateIncrements(Ui, ui); }
void FETiedFluidSurface::PrepStep() { m_lc.PrepStep(); }

//-----------------------------------------------------------------------------
// FETiedFluidInterface
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Define sliding interface parameters
BEGIN_FECORE_CLASS(FETiedFluidInterface, FEContactInterface)
    ADD_PARAMETER(m_laugon   , "laugon")->setLongName("Enforcement method")->setEnums("PENALTY\0AUGLAG\0");
    ADD_PARAMETER(m_tol      , "tolerance"          );
    ADD_PARAMETER(m_eps      , "penalty"            );
    ADD_PARAMETER(m_btwopass , "two_pass"           );
    ADD_PARAMETER(m_stol     , "search_tol"         );
    ADD_PARAMETER(m_srad     , "search_radius"      )->setUnits(UNIT_LENGTH);
    ADD_PARAMETER(m_naugmin  , "minaug"             );
    ADD_PARAMETER(m_naugmax  , "maxaug"             );
END_FECORE_CLASS();

FETiedFluidInterface::FETiedFluidInterface(FEModel* pfem) : FEContactInterface(pfem), m_ss(pfem), m_ms(pfem), m_dofWE(pfem)
{
    static int count = 1;
    SetID(count++);
    
    // initial values
    m_tol = 0.1;
    m_eps = 1;
    m_stol = 0.01;
    m_srad = 1.0;
    m_btwopass = true;
    m_bfreedofs = true;
    
    m_naugmin = 0;
    m_naugmax = 10;
    
    // set parents
    m_ss.SetContactInterface(this);
    m_ms.SetContactInterface(this);

    m_ss.SetSibling(&m_ms);
    m_ms.SetSibling(&m_ss);
}

//-----------------------------------------------------------------------------

FETiedFluidInterface::~FETiedFluidInterface()
{
}

//-----------------------------------------------------------------------------
bool FETiedFluidInterface::Init()
{
    // initialize surface data
    if (m_ss.Init() == false) return false;
    if (m_ms.Init() == false) return false;
    
    // get the DOFS
    m_dofWE.AddVariable(FEBioFluid::GetVariableName(FEBioFluid::RELATIVE_FLUID_VELOCITY));
    m_dofWE.AddVariable(FEBioFluid::GetVariableName(FEBioFluid::FLUID_DILATATION));
    
    return true;
}

//-----------------------------------------------------------------------------
void FETiedFluidInterface::Activate()
{
    // don't forget to call the base class
    FEContactInterface::Activate();
    
    // project surfaces onto each other (node-to-surface projection)
    InitialNodalProjection(m_ss, m_ms);
    InitialNodalProjection(m_ms, m_ss);
    
    m_ss.Activate();
    if (m_btwopass) m_ms.Activate();

}

//-----------------------------------------------------------------------------
// Perform initial projection between tied surfaces in reference configuration
void FETiedFluidInterface::InitialNodalProjection(FETiedFluidSurface& ss, FETiedFluidSurface& ms)
{
    // get primary surface
    FESurface& psurf = *ss.GetSurface();
    psurf.UpdateNodeNormals();
    
    // get secondary surface
    FESurface& ssurf = *ms.GetSurface();
    ssurf.UpdateNodeNormals();
    
    FEMesh& mesh = GetMesh();
    FESurfaceElement* pme;
    vec3d r, nu;
    double rs[2];
    
    // initialize projection data
    FENormalProjection np(ssurf);
    np.SetTolerance(m_stol);
    np.SetSearchRadius(m_srad);
    np.Init();
    
    
    // loop over all nodes
    int n = 0;
    for (int i=0; i<psurf.Nodes(); ++i)
    {
        FENode& node = psurf.Node(i);
        r = node.m_rt;
            
        // calculate the normal at this integration point
        nu = psurf.NodeNormal(i);
        
        // find the intersection point with the secondary surface
        pme = np.Project2(r, nu, rs);
        
        ss.m_pme[i] = pme;
        ss.m_rs[i][0] = rs[0];
        ss.m_rs[i][1] = rs[1];
    }
}

//-----------------------------------------------------------------------------
void FETiedFluidInterface::Update(const std::vector<double>& Ui, const std::vector<double>& ui)
{
    m_ss.Update(Ui, ui);
    if (m_btwopass) m_ms.Update(Ui, ui);
}

//-----------------------------------------------------------------------------
void FETiedFluidInterface::UpdateIncrements(std::vector<double>& Ui, const std::vector<double>& ui)
{
    m_ss.UpdateIncrements(Ui, ui);
    if (m_btwopass) m_ms.UpdateIncrements(Ui, ui);
}

//-----------------------------------------------------------------------------
void FETiedFluidInterface::LoadVector(FEGlobalVector& R, const FETimeInfo& tp)
{
    m_ss.LoadVector(R, tp);
    if (m_btwopass) m_ms.LoadVector(R, tp);
}

//-----------------------------------------------------------------------------
void FETiedFluidInterface::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
{
    m_ss.StiffnessMatrix(LS, tp);
    if (m_btwopass) m_ms.StiffnessMatrix(LS, tp);
}

//-----------------------------------------------------------------------------
void FETiedFluidInterface::BuildMatrixProfile(FEGlobalMatrix& G)
{
    m_ss.BuildMatrixProfile(G);
    if (m_btwopass) m_ms.BuildMatrixProfile(G);
}

//-----------------------------------------------------------------------------
int FETiedFluidInterface::InitEquations(int neq)
{
    m_ss.InitEquations(neq);
    if (m_btwopass) m_ms.InitEquations(neq);
}

//-----------------------------------------------------------------------------
bool FETiedFluidInterface::Augment(int naug, const FETimeInfo& tp)
{
    if (!m_ss.Augment(naug, tp)) return false;
    if (m_btwopass) {
        if (!m_ms.Augment(naug, tp)) return false;
    }
    return true;
}

//-----------------------------------------------------------------------------
void FETiedFluidInterface::PrepStep()
{
    m_ss.PrepStep();
    if (m_btwopass) m_ms.PrepStep();
}

//-----------------------------------------------------------------------------
void FETiedFluidInterface::Serialize(DumpStream &ar)
{
    // store contact data
    FEContactInterface::Serialize(ar);
    
    // store contact surface data
    m_ss.Serialize(ar);
    m_ms.Serialize(ar);
    
    if (ar.IsShallow()) return;
    ar & m_dofWE;
}
