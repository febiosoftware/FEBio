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
#include "FEFluidFSI.h"

//-----------------------------------------------------------------------------
// FETiedFluidSurface
//-----------------------------------------------------------------------------

FETiedFluidSurface::FETiedFluidSurface(FEModel* pfem) : FESurface(pfem)
{
}

//-----------------------------------------------------------------------------
bool FETiedFluidSurface::Init()
{
    // initialize surface data first
    if (FESurface::Init() == false) return false;
    
    FESurface* surf = GetSurface();
    m_pme.assign(surf->Nodes(), nullptr);
    double rs[2] = {0,0};
    m_rs.assign(surf->Nodes(), std::vector<double>(2,0));
    m_gap.assign(surf->Nodes(), 0.0);
    m_tag.assign(surf->Nodes(), -1);

    return true;
}

//-----------------------------------------------------------------------------
// FETiedFluidInterface
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Define sliding interface parameters
BEGIN_FECORE_CLASS(FETiedFluidInterface, FEContactInterface)
    ADD_PARAMETER(m_lcv.m_laugon   , "velocity_laugon")->setLongName("Enforcement method")->setEnums("PENALTY\0AUGLAG\0LAGMULT\0");
    ADD_PARAMETER(m_lcv.m_tol      , "velocity_tol"       );
    ADD_PARAMETER(m_lcv.m_eps      , "velocity_penalty"   );
    ADD_PARAMETER(m_lcv.m_naugmin  , "velocity_minaug"    );
    ADD_PARAMETER(m_lcv.m_naugmax  , "velocity_maxaug"    );

    ADD_PARAMETER(m_lcp.m_laugon   , "pressure_laugon")->setLongName("Enforcement method")->setEnums("PENALTY\0AUGLAG\0LAGMULT\0");
    ADD_PARAMETER(m_lcp.m_tol      , "pressure_tol"       );
    ADD_PARAMETER(m_lcp.m_eps      , "pressurepenalty"   );
    ADD_PARAMETER(m_lcp.m_naugmin  , "pressure_minaug"    );
    ADD_PARAMETER(m_lcp.m_naugmax  , "pressure_maxaug"    );

    ADD_PARAMETER(m_btwopass , "two_pass"           );
    ADD_PARAMETER(m_stol     , "search_tol"         );
    ADD_PARAMETER(m_srad     , "search_radius"      )->setUnits(UNIT_LENGTH);
    ADD_PARAMETER(m_bfreedofs, "free_DOFs")->setLongName("Free fixed interface DOFs");
END_FECORE_CLASS();

FETiedFluidInterface::FETiedFluidInterface(FEModel* pfem) : FEContactInterface(pfem), m_ss(pfem), m_ms(pfem), m_dofWE(pfem), m_lcv(pfem), m_lcp(pfem)
{
    static int count = 1;
    SetID(count++);
    
    // initial values
    m_lcv.m_laugon = m_lcp.m_laugon = 0;
    m_lcv.m_tol = m_lcp.m_tol = 0.1;
    m_lcv.m_eps = m_lcp.m_eps = 1;
    m_lcv.m_tol = 0.01;
    m_lcv.m_naugmin = m_lcp.m_naugmin = 0;
    m_lcv.m_naugmax = m_lcp.m_naugmax = 10;
    
    m_stol = 0.01;
    m_srad = 1.0;
    m_btwopass = true;
    m_bfreedofs = true;
    
    m_binit = false;
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
    
//    m_laugon = max(m_lcv.m_laugon, m_lcp.m_laugon);
    
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
    
    if (m_binit == false)
    {
        int npass = (m_btwopass?2:1);
        for (int np=0; np<npass; ++np)
        {
            FETiedFluidSurface& ps = (np == 0? m_ss : m_ms);
            FETiedFluidSurface& ss = (np == 0? m_ms : m_ss);

            // don't forget to call base class
            
            FEModel& fem = *GetFEModel();
            FEMesh& mesh = fem.GetMesh();
            FESurface* surf = ss.GetSurface();
            
            // create linear constraints
            for (int i=0; i<ps.Nodes(); ++i) {
                FENode& node = ps.Node(i);
                if (ps.m_pme[i] == nullptr) continue;
                if ((node.HasFlags(FENode::EXCLUDE) == false) && (node.m_rid == -1)) {
                    
                    // get shape functions on sibling surface at projection point
                    double Na[FEElement::MAX_NODES];
                    ps.m_pme[i]->shape_fnc(Na, ps.m_rs[i][0], ps.m_rs[i][1]);
                    
                    // free the fixed DOFs inside of the tied interface, if requested
                    if (m_bfreedofs) {
                        for (int j=0; j<m_dofWE.Size(); ++j) {
                            // check to see if the DOF is not open
                            int dof = m_dofWE[j];
                            if (node.get_bc(dof) == DOF_FIXED) {
                                // if not open, check to see if one of the projected face nodal DOFs is not open
                                bool reset = true;
                                for (int k=0; k<ps.m_pme[i]->Nodes(); ++k) {
                                    FENode& n2 = surf->Node(ps.m_pme[i]->m_lnode[k]);
                                    int indx = ps.m_pme[i]->m_lnode[k];
                                    if ((ss.m_tag[indx] == -1) && (n2.get_bc(dof) == DOF_FIXED) && (fabs(Na[k]) > m_stol)) {
                                        ss.m_tag[indx] = 1;
                                        reset = false;
                                    }
                                }
                                // only reset if there are no significant constraints on any node of the projected surface
                                if (reset) node.set_bc(dof, DOF_OPEN);
                            }
                            else if (node.get_bc(dof) == DOF_PRESCRIBED) {
                                feLogError("Prescribed degrees of freedom should not be assigned within tied-fluid interfaces!");
                                exit(-1);
                            }
                        }
                    }
                    
                    // extract the fluid bulk modulus for the solid element underneath this face
                    FEElement& el = *(ps.m_pme[i]->m_elem[0].pe);
                    FEMaterial* pm = fem.GetMaterial(el.GetMatID());
                    double K = 0;
                    FEFluidMaterial* fm = dynamic_cast<FEFluidMaterial*>(pm);
                    FEFluidFSI* fsim = dynamic_cast<FEFluidFSI*>(pm);
                    if (fm) K = fm->BulkModulus(*el.GetMaterialPoint(0));
                    else if (fsim) K = fsim->Fluid()->BulkModulus(*el.GetMaterialPoint(0));
                    
                    //  constrain components of velocities on primary and secondary surfaces
                    FEAugLagLinearConstraint *pLCvx = fecore_alloc(FEAugLagLinearConstraint, GetFEModel());
                    FEAugLagLinearConstraint *pLCvy = fecore_alloc(FEAugLagLinearConstraint, GetFEModel());
                    FEAugLagLinearConstraint *pLCvz = fecore_alloc(FEAugLagLinearConstraint, GetFEModel());
                    FEAugLagLinearConstraint *pLCef = fecore_alloc(FEAugLagLinearConstraint, GetFEModel());
                    pLCvx->AddDOF(node.GetID(), m_dofWE[0], 1);
                    pLCvy->AddDOF(node.GetID(), m_dofWE[1], 1);
                    pLCvz->AddDOF(node.GetID(), m_dofWE[2], 1);
                    pLCef->AddDOF(node.GetID(), m_dofWE[3], 1*K);
                    for (int j=0; j<ps.m_pme[i]->Nodes(); ++j) {
                        FENode& n2 = ss.Node(ps.m_pme[i]->m_lnode[j]);
                        pLCvx->AddDOF(n2.GetID(), m_dofWE[0], -Na[j]);
                        pLCvy->AddDOF(n2.GetID(), m_dofWE[1], -Na[j]);
                        pLCvz->AddDOF(n2.GetID(), m_dofWE[2], -Na[j]);
                        pLCef->AddDOF(n2.GetID(), m_dofWE[3], -Na[j]*K);
                    }
                    m_lcv.add(pLCvx);
                    m_lcv.add(pLCvy);
                    m_lcv.add(pLCvz);
                    m_lcp.add(pLCef);
                }
            }
        }
        m_lcv.Init();
        m_lcv.Activate();
        m_lcp.Init();
        m_lcp.Activate();
        m_binit = true;
    }
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
    double rs[2] = {0};
    
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
        if (pme != nullptr) {
            // evaluate the location of the intersection point in 3D
            vec3d r2[FEElement::MAX_NODES];
            for (int j=0; j<pme->Nodes(); ++j)
                r2[j]= mesh.Node(pme->m_node[j]).m_rt;
            vec3d q = pme->eval(r2,rs[0], rs[1]);
            ss.m_rs[i][0] = rs[0];
            ss.m_rs[i][1] = rs[1];
            ss.m_gap[i] = (q - r).norm();
        }
    }
}

//-----------------------------------------------------------------------------
void FETiedFluidInterface::Update(std::vector<double>& Ui, std::vector<double>& ui)
{
    m_lcv.Update(Ui, ui);
    m_lcp.Update(Ui, ui);
}

//-----------------------------------------------------------------------------
void FETiedFluidInterface::UpdateIncrements(std::vector<double>& Ui, const std::vector<double>& ui)
{
    m_lcv.UpdateIncrements(Ui, ui);
    m_lcp.UpdateIncrements(Ui, ui);
}

//-----------------------------------------------------------------------------
void FETiedFluidInterface::LoadVector(FEGlobalVector& R, const FETimeInfo& tp)
{
    m_lcv.LoadVector(R, tp);
    m_lcp.LoadVector(R, tp);
}

//-----------------------------------------------------------------------------
void FETiedFluidInterface::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
{
    m_lcv.StiffnessMatrix(LS, tp);
    m_lcp.StiffnessMatrix(LS, tp);
}

//-----------------------------------------------------------------------------
void FETiedFluidInterface::BuildMatrixProfile(FEGlobalMatrix& G)
{
    m_lcv.BuildMatrixProfile(G);
    m_lcp.BuildMatrixProfile(G);
}

//-----------------------------------------------------------------------------
int FETiedFluidInterface::InitEquations(int neq)
{
    int n = m_lcv.InitEquations(neq);
    n += m_lcp.InitEquations(neq);
	return n;
}

//-----------------------------------------------------------------------------
bool FETiedFluidInterface::Augment(int naug, const FETimeInfo& tp)
{
    if (!m_lcv.Augment(naug, tp) || !m_lcp.Augment(naug, tp)) return false;
    return true;
}

//-----------------------------------------------------------------------------
void FETiedFluidInterface::PrepStep()
{
    m_lcv.PrepStep();
    m_lcp.PrepStep();
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
