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
#include "FEThermoFluidTemperatureBC.h"
#include "FEBioThermoFluid.h"
#include "FEThermoFluid.h"
#include <FECore/FEModel.h>

//-----------------------------------------------------------------------------
//! constructor
FEThermoFluidTemperatureBC::FEThermoFluidTemperatureBC(FEModel* pfem) : FEPrescribedSurface(pfem), m_dofW(pfem)
{
    m_dofT = -1;
    m_surf = nullptr;
}

//-----------------------------------------------------------------------------
//! initialize
bool FEThermoFluidTemperatureBC::Init()
{
    FEModel& fem = *GetFEModel();
    
    m_dofW.Clear();
    m_dofW.AddVariable(FEBioThermoFluid::GetVariableName(FEBioThermoFluid::RELATIVE_FLUID_VELOCITY));
    m_dofT = fem.GetDOFIndex(FEBioThermoFluid::GetVariableName(FEBioThermoFluid::TEMPERATURE), 0);
    SetDOFList(m_dofW);
    SetDOFList(m_dofT);

    if (FEPrescribedSurface::Init() == false) return false;

    m_surf = GetSurface();
    m_T.assign(m_surf->Nodes(), 0.0);
    m_backflow.assign(m_surf->Nodes(), false);
    
    m_nnlist.Create(fem.GetMesh());
    
    return FEPrescribedSurface::Init();
}

//-----------------------------------------------------------------------------
//! Evaluate and prescribe the resistance pressure
void FEThermoFluidTemperatureBC::Update()
{
    // determine backflow conditions
    MarkBackFlow();
    
    int N = m_surf->Nodes();
    std::vector<vector<double>> efNodes(N, vector<double>());
    std::vector<vector<double>> vn(N, vector<double>());
    
    for (int i=0; i<m_surf->Nodes(); ++i)
    {
        FENode& node = m_surf->Node(i);
        
        if (m_backflow[i]) m_T[i] = node.get_prev(m_dofT);
        else {
            int nid = node.GetID()-1; //0 based
            int val = m_nnlist.Valence(nid);
            int* nlist = m_nnlist.NodeList(nid);
            
            vector<int> connectnidarray;
            int connectnid=0;
            //check which connecting nodes are not on surface
            for (int j = 0; j<val; ++j)
            {
                int cnid = *(nlist+j);
                bool isSurf = false;
                for (int k=0; k<m_surf->Nodes(); ++k)
                {
                    if (cnid==m_surf->Node(k).GetID()-1)
                        isSurf = true;
                }
                if (!isSurf)
                {
                    connectnidarray.push_back(cnid);
                }
            }
            //find closest connecting node
            if (connectnidarray.size()>0)
            {
                int cnodeIndex = 0;
                FENode* cnodeArray = GetFEModel()->GetMesh().FindNodeFromID(connectnidarray[cnodeIndex]+1);
                double smallDist = sqrt(pow((node.m_rt.x-cnodeArray->m_rt.x),2)+pow((node.m_rt.y-cnodeArray->m_rt.y),2)+pow((node.m_rt.z-cnodeArray->m_rt.z),2));
                for (int j = 1; j<connectnidarray.size(); ++j)
                {
                    FENode* tempNode = GetFEModel()->GetMesh().FindNodeFromID(connectnidarray[j]+1);
                    double temp = sqrt(pow((node.m_rt.x-tempNode->m_rt.x),2)+pow((node.m_rt.y-tempNode->m_rt.y),2)+pow((node.m_rt.z-tempNode->m_rt.z),2));
                    if(temp<smallDist)
                    {
                        cnodeIndex = j;
                        cnodeArray = tempNode;
                        smallDist = temp;
                    }
                }
                connectnid = connectnidarray[cnodeIndex];
                FENode* cnode = GetFEModel()->GetMesh().FindNodeFromID(connectnid+1);
                m_T[i] = cnode->get(m_dofT);
            }
            else
                m_T[i] = node.get_prev(m_dofT);
        }
    }
    FEPrescribedSurface::Update();
}

//-----------------------------------------------------------------------------
//! evaluate the flow rate across this surface
void FEThermoFluidTemperatureBC::MarkBackFlow()
{
    // Calculate normal flow velocity on each face to determine
    // backflow condition
    vec3d rt[FEElement::MAX_NODES];
    vec3d vt[FEElement::MAX_NODES];
    
    const FETimeInfo& tp = GetTimeInfo();
    
    for (int iel=0; iel<m_surf->Elements(); ++iel)
    {
        FESurfaceElement& el = m_surf->Element(iel);
        
        // nr integration points
        int nint = el.GaussPoints();
        
        // nr of element nodes
        int neln = el.Nodes();
        
        // nodal coordinates
        for (int i=0; i<neln; ++i) {
            FENode& node = m_surf->GetMesh()->Node(el.m_node[i]);
            rt[i] = node.m_rt*tp.alpha + node.m_rp*(1-tp.alpha);
            vt[i] = node.get_vec3d(m_dofW[0], m_dofW[1], m_dofW[2])*tp.alpha + node.get_vec3d_prev(m_dofW[0], m_dofW[1], m_dofW[2])*(1-tp.alpha);
        }
        
        double* Nr, *Ns;
        double* N;
        double* w  = el.GaussWeights();
        
        vec3d dxr, dxs, v;
        double vn = 0;
        
        // repeat over integration points
        for (int n=0; n<nint; ++n)
        {
            N  = el.H(n);
            Nr = el.Gr(n);
            Ns = el.Gs(n);
            
            // calculate the velocity and tangent vectors at integration point
            dxr = dxs = v = vec3d(0,0,0);
            for (int i=0; i<neln; ++i)
            {
                v += vt[i]*N[i];
                dxr += rt[i]*Nr[i];
                dxs += rt[i]*Ns[i];
            }
            
            vec3d normal = dxr ^ dxs;
            normal.unit();
            vn += (normal*v)*w[n];
        }
        
        if (vn < 0)
            for (int i=0; i<neln; ++i) m_backflow[el.m_lnode[i]] = true;
        else
            for (int i=0; i<neln; ++i) m_backflow[el.m_lnode[i]] = false;

    }
}

//-----------------------------------------------------------------------------
void FEThermoFluidTemperatureBC::GetNodalValues(int nodelid, std::vector<double>& val)
{
    val[0] = m_T[nodelid];
}

//-----------------------------------------------------------------------------
// copy data from another class
void FEThermoFluidTemperatureBC::CopyFrom(FEBoundaryCondition* pbc)
{
    // TODO: implement this
    assert(false);
}

//-----------------------------------------------------------------------------
//! serialization
void FEThermoFluidTemperatureBC::Serialize(DumpStream& ar)
{
    FEPrescribedSurface::Serialize(ar);
    ar & m_T;
    if (ar.IsShallow()) return;
    ar & m_dofT;
    ar & m_surf;
}
