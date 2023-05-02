//
//  FEFluidSolutesPressure.cpp
//  FEBioFluid
//
//  Created by Jay Shim on 12/10/20.
//  Copyright © 2020 febio.org. All rights reserved.
//

#include "stdafx.h"
#include "FEFluidSolutesPressure.h"
#include "FEBioFluidSolutes.h"
#include <FECore/FEModel.h>

//=============================================================================
BEGIN_FECORE_CLASS(FEFluidSolutesPressure, FESurfaceLoad)
ADD_PARAMETER(m_p, "pressure");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FEFluidSolutesPressure::FEFluidSolutesPressure(FEModel* pfem) : FESurfaceLoad(pfem)
{
    m_pfs = nullptr;
    m_p = 0;
    

    m_dofEF = pfem->GetDOFIndex(FEBioFluidSolutes::GetVariableName(FEBioFluidSolutes::FLUID_DILATATION), 0);
    m_dofC = pfem->GetDOFIndex(FEBioFluidSolutes::GetVariableName(FEBioFluidSolutes::FLUID_CONCENTRATION), 0);
    
    m_dof.Clear();
    m_dof.AddDof(m_dofEF);
    
}

//-----------------------------------------------------------------------------
//! initialize
bool FEFluidSolutesPressure::Init()
{
    if (FESurfaceLoad::Init() == false) return false;
    
    // get fluid from first surface element
    // assuming the entire surface bounds the same fluid
    FESurfaceElement& el = m_psurf->Element(0);
    FEElement* pe = el.m_elem[0];
    if (pe == nullptr) return false;
    
    // get the material
    FEMaterial* pm = GetFEModel()->GetMaterial(pe->GetMatID());
    m_pfs = dynamic_cast<FEFluidSolutes*>(pm);
    if (m_pfs == nullptr) return false;
    
    return true;
}

//-----------------------------------------------------------------------------
//! Activate the degrees of freedom for this BC
void FEFluidSolutesPressure::Activate()
{
    FESurface* ps = &GetSurface();
    
    for (int i=0; i<ps->Nodes(); ++i)
    {
        FENode& node = ps->Node(i);
        // mark node as having prescribed DOF
        node.set_bc(m_dofEF, DOF_PRESCRIBED);
    }
    
    FESurfaceLoad::Activate();
}

//-----------------------------------------------------------------------------
//! Evaluate and prescribe the resistance pressure
void FEFluidSolutesPressure::Update()
{
    // prescribe this dilatation at the nodes
    FESurface* ps = &GetSurface();
    int nsol = 0;
    double T = 0;
    double R = 0;
    if (m_pfs){
        nsol = m_pfs->Solutes();
        T = m_pfs->m_Tabs;
        R = m_pfs->m_Rgas;
    }
    
    std::map<int,vector<double>> oscNodes;
    std::map<int,vector<double>> caNodes;
    
    //initialize maps based on surface nodal IDs
    for (int i=0; i<ps->Nodes(); ++i)
    {
        oscNodes.insert(pair<int,vector<double> >(ps->Node(i).GetID(), vector<double>()));
        caNodes.insert(pair<int,vector<double> >(ps->Node(i).GetID(), vector<double>()));
    }
    
    //Project sum of all ca and osc values from int points to nodes on surface
    //All values put into map, including duplicates
    for (int i=0; i<ps->Elements(); ++i)
    {
        FESurfaceElement& el = ps->Element(i);
        FEElement* e = el.m_elem[0];
        FESolidElement* se = dynamic_cast<FESolidElement*>(e);
        if (se) {
            double osci[FEElement::MAX_INTPOINTS];
            double osco[FEElement::MAX_NODES];
            double cai[FEElement::MAX_INTPOINTS];
            double cao[FEElement::MAX_NODES];
            for (int j=0; j<se->GaussPoints(); ++j) {
                FEMaterialPoint* pt = se->GetMaterialPoint(j);
                FEFluidSolutesMaterialPoint* fsp = pt->ExtractData<FEFluidSolutesMaterialPoint>();
                if (fsp)
                {
                    if (m_pfs)
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
            for (int j=0; j<el.Nodes(); ++j)
            {
                oscNodes[el.m_node[j]+1].push_back(osco[se->FindNode(el.m_node[j])]);
                caNodes[el.m_node[j]+1].push_back(cao[se->FindNode(el.m_node[j])]);
            }
        }
        //If no solid element, insert all 0s
        else{
            for (int j=0; j<el.Nodes(); ++j)
            {
                oscNodes[el.m_node[j]+1].push_back(0);
                caNodes[el.m_node[j]+1].push_back(0);
            }
        }
    }
    //For each node average the nodal ca and osc and then calculate ef based on desired p
    for (int i=0; i<ps->Nodes(); ++i)
    {
        if (ps->Node(i).m_ID[m_dofEF] < -1)
        {
            FENode& node = ps->Node(i);
            
            //get osmotic component of pressure
            double ca = 0;
            double osc = 0;
            for (int j = 0; j < caNodes[node.GetID()].size(); ++j)
            {
                ca += caNodes[node.GetID()][j];
                osc += oscNodes[node.GetID()][j];
            }
            ca /= caNodes[node.GetID()].size();
            osc /= caNodes[node.GetID()].size();
            
            double pc = R*T*osc*ca;
            
            //get correct ef for desired pressure
            double e = 0;
            bool good = false;
            if (m_pfs)
                good = m_pfs->Fluid()->Dilatation(0, m_p - pc, e);
            assert(good);
            
            // set node as having prescribed DOF
            node.set(m_dofEF, e);
        }
    }
    
    /*
    for (int i=0; i<ps->Nodes(); ++i)
    {
        if (ps->Node(i).m_ID[m_dofEF] < -1)
        {
            FENode& node = ps->Node(i);
            
            //get osmotic component of pressure
            double c = 0;
            for (int isol = 0; isol < nsol; ++isol)
            {
                c += node.get(m_dofC + isol);
            }
            
            //get correct ef for desired pressure
            double e = -1/K*(m_p - R*T*m_osc*c);
            
            // set node as having prescribed DOF
            node.set(m_dofEF, e);
        }
    }
    */
    GetFEModel()->SetMeshUpdateFlag(true);
}

//-----------------------------------------------------------------------------
//! calculate residual
void FEFluidSolutesPressure::LoadVector(FEGlobalVector& R)
{
}

//-----------------------------------------------------------------------------
//! serialization
void FEFluidSolutesPressure::Serialize(DumpStream& ar)
{
    FESurfaceLoad::Serialize(ar);
    if (ar.IsShallow()) return;
    ar & m_pfs;
    ar & m_dofC & m_dofEF;
}
