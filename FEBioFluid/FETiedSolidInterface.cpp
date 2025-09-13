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
#include "FETiedSolidInterface.h"
#include <FECore/FENormalProjection.h>
#include <FECore/FEClosestPointProjection.h>
#include <FECore/log.h>
#include <FECore/DumpStream.h>
#include <FECore/FEGlobalMatrix.h>
#include <FECore/FELinearSystem.h>
#include <FECore/FEModel.h>
#include "FEBioFSI.h"
#include "FEFluidFSI.h"

//-----------------------------------------------------------------------------
// FETiedFluidFSISurface
//-----------------------------------------------------------------------------

FETiedFluidFSISurface::FETiedFluidFSISurface(FEModel* pfem) : FESurface(pfem)
{
}

//-----------------------------------------------------------------------------
bool FETiedFluidFSISurface::Init()
{
    // initialize surface data first
    if (FESurface::Init() == false) return false;
    
    FESurface* surf = GetSurface();
    m_pme.assign(surf->Nodes(), nullptr);
    double rs[2] = {0,0};
    m_rs.assign(surf->Nodes(), std::vector<double>(2,0));
    
    m_tag.assign(surf->Nodes(),-1);

    return true;
}
//-----------------------------------------------------------------------------
// FETiedSolidInterface
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Define sliding interface parameters
BEGIN_FECORE_CLASS(FETiedSolidInterface, FEContactInterface)
    ADD_PARAMETER(m_lc.m_laugon   , "laugon")->setLongName("Enforcement method")->setEnums("PENALTY\0AUGLAG\0LAGMULT\0");
    ADD_PARAMETER(m_lc.m_tol      , "tolerance"          );
    ADD_PARAMETER(m_lc.m_eps      ,"displacement_penalty");
    ADD_PARAMETER(m_lc.m_naugmin  , "minaug"             );
    ADD_PARAMETER(m_lc.m_naugmax  , "maxaug"             );

    ADD_PARAMETER(m_btwopass , "two_pass"           );
    ADD_PARAMETER(m_stol     , "search_tol"         );
    ADD_PARAMETER(m_srad     , "search_radius"      )->setUnits(UNIT_LENGTH);
    ADD_PARAMETER(m_bfreedofs, "free_DOFs")->setLongName("Free fixed interface DOFs");
END_FECORE_CLASS();

FETiedSolidInterface::FETiedSolidInterface(FEModel* pfem) : FEContactInterface(pfem), m_ss(pfem), m_ms(pfem), m_dofU(pfem), m_lc(pfem)
{
    static int count = 1;
    SetID(count++);
    
    // initial values
    m_lc.m_laugon = 0;
    m_lc.m_tol = 0.1;
    m_lc.m_eps = 1;
    m_lc.m_tol = 0.01;
    m_lc.m_naugmin = 0;
    m_lc.m_naugmax = 10;
    
    m_stol = 0.01;
    m_srad = 1.0;
    m_btwopass = true;
    m_bfreedofs = false;
    
    m_binit = false;
}

//-----------------------------------------------------------------------------

FETiedSolidInterface::~FETiedSolidInterface()
{
}

//-----------------------------------------------------------------------------
bool FETiedSolidInterface::Init()
{
    // initialize surface data
    if (m_ss.Init() == false) return false;
    if (m_ms.Init() == false) return false;

    m_laugon = m_lc.m_laugon;

    // get the DOFS
    m_dofU.AddVariable(FEBioFSI::GetVariableName(FEBioFSI::DISPLACEMENT));
    
    return true;
}

//-----------------------------------------------------------------------------
void FETiedSolidInterface::Activate()
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
            FETiedFluidFSISurface& ps = (np == 0? m_ss : m_ms);
            FETiedFluidFSISurface& ss = (np == 0? m_ms : m_ss);

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
                        for (int j=0; j<m_dofU.Size(); ++j) {
                            // check to see if the DOF is not open
                            int dof = m_dofU[j];
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
                    FEAugLagLinearConstraint *pLCux = fecore_alloc(FEAugLagLinearConstraint, GetFEModel());
                    FEAugLagLinearConstraint *pLCuy = fecore_alloc(FEAugLagLinearConstraint, GetFEModel());
                    FEAugLagLinearConstraint *pLCuz = fecore_alloc(FEAugLagLinearConstraint, GetFEModel());
                    pLCux->AddDOF(node.GetID(), m_dofU[0], 1);
                    pLCuy->AddDOF(node.GetID(), m_dofU[1], 1);
                    pLCuz->AddDOF(node.GetID(), m_dofU[2], 1);
                    for (int j=0; j<ps.m_pme[i]->Nodes(); ++j) {
                        FENode& n2 = ss.Node(ps.m_pme[i]->m_lnode[j]);
                        pLCux->AddDOF(n2.GetID(), m_dofU[0], -Na[j]);
                        pLCuy->AddDOF(n2.GetID(), m_dofU[1], -Na[j]);
                        pLCuz->AddDOF(n2.GetID(), m_dofU[2], -Na[j]);
                    }
                    m_lc.add(pLCux);
                    m_lc.add(pLCuy);
                    m_lc.add(pLCuz);
                }
            }
            m_lc.Init();
            m_lc.Activate();
            m_binit = true;
        }
    }
}

//-----------------------------------------------------------------------------
// Perform initial projection between tied surfaces in reference configuration
void FETiedSolidInterface::InitialNodalProjection(FETiedFluidFSISurface& ss, FETiedFluidFSISurface& ms)
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
void FETiedSolidInterface::Update(const std::vector<double>& Ui, const std::vector<double>& ui)
{
    m_lc.Update(Ui, ui);
}

//-----------------------------------------------------------------------------
void FETiedSolidInterface::UpdateIncrements(std::vector<double>& Ui, const std::vector<double>& ui)
{
    m_lc.UpdateIncrements(Ui, ui);
}

//-----------------------------------------------------------------------------
void FETiedSolidInterface::LoadVector(FEGlobalVector& R, const FETimeInfo& tp)
{
    m_lc.LoadVector(R, tp);
}

//-----------------------------------------------------------------------------
void FETiedSolidInterface::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
{
    m_lc.StiffnessMatrix(LS, tp);
}

//-----------------------------------------------------------------------------
void FETiedSolidInterface::BuildMatrixProfile(FEGlobalMatrix& G)
{
    m_lc.BuildMatrixProfile(G);
}

//-----------------------------------------------------------------------------
int FETiedSolidInterface::InitEquations(int neq)
{
	return m_lc.InitEquations(neq);
}

//-----------------------------------------------------------------------------
bool FETiedSolidInterface::Augment(int naug, const FETimeInfo& tp)
{
    if (!m_lc.Augment(naug, tp)) return false;
    return true;
}

//-----------------------------------------------------------------------------
void FETiedSolidInterface::PrepStep()
{
    m_lc.PrepStep();
}

//-----------------------------------------------------------------------------
void FETiedSolidInterface::Serialize(DumpStream &ar)
{
    // store contact data
    FEContactInterface::Serialize(ar);
    
    // store contact surface data
    m_ss.Serialize(ar);
    m_ms.Serialize(ar);
    
    if (ar.IsShallow()) return;
    ar & m_dofU;
}
