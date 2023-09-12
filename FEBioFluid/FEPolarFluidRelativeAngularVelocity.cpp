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
#include "FEPolarFluidRelativeAngularVelocity.h"
#include "FEPolarFluid.h"
#include <FECore/FEAnalysis.h>
#include <FECore/FEModel.h>
#include <FECore/FESurface.h>
#include "FEBioPolarFluid.h"

//=============================================================================
BEGIN_FECORE_CLASS(FEPolarFluidRelativeAngularVelocity, FEPrescribedSurface)
    ADD_PARAMETER(m_h, "h")->setUnits(UNIT_ANGULAR_VELOCITY)->setLongName("relative angular velocity");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FEPolarFluidRelativeAngularVelocity::FEPolarFluidRelativeAngularVelocity(FEModel* pfem) : FEPrescribedSurface(pfem), m_dofG(pfem)
{
    m_h = vec3d(0,0,0);
    if (pfem)
    {
        m_dofG.AddVariable(FEBioPolarFluid::GetVariableName(FEBioPolarFluid::FLUID_ANGULAR_VELOCITY));
    }
    m_psurf = nullptr;
}

//-----------------------------------------------------------------------------
//! initialize
bool FEPolarFluidRelativeAngularVelocity::Init()
{
    SetDOFList(m_dofG);

    if (FEPrescribedSurface::Init() == false) return false;
    
    m_psurf = FESurfaceBC::GetSurface();

    m_g.assign(GetSurface()->Nodes(), vec3d(0,0,0));

    return true;
}

//-----------------------------------------------------------------------------
//! Evaluate and prescribe the resistance pressure
void FEPolarFluidRelativeAngularVelocity::Update()
{
    // prescribe this dilatation at the nodes
    FESurface* ps = GetSurface();
    FEModel* fem = GetFEModel();
    FEMesh& mesh = fem->GetMesh();
    FETimeInfo& tp = fem->GetTime();

    int N = ps->Nodes();
    std::vector< vector<vec3d>> gNodes(N);

    //Project sum of all ca and osc values from int points to nodes on surface
    //All values put into map, including duplicates
    for (int i=0; i<ps->Elements(); ++i)
    {
        FESurfaceElement& el = ps->Element(i);
        // evaluate average prescribed relative angular velocity on this face
        vec3d h = vec3d(0,0,0);
        for (int j=0; j<el.GaussPoints(); ++j) {
            FEMaterialPoint* pt = el.GetMaterialPoint(j);
            h += m_h(*pt);
        }
        h /= el.GaussPoints();
        // get surface underlying material
        FEElement* e = el.m_elem[0];
        FESolidElement* se = dynamic_cast<FESolidElement*>(e);
        if (se) {
            FEMaterial* pm = GetFEModel()->GetMaterial(e->GetMatID());
            FEPolarFluid* pofl = pm->ExtractProperty<FEPolarFluid>();
            if (pofl) {
                // evalauate average regional angular velocity w
                vec3d w(0,0,0);
                for (int j=0; j<se->GaussPoints(); ++j) {
                    FEMaterialPoint* pt = se->GetMaterialPoint(j);
                    FEPolarFluidMaterialPoint* ppt = pt->ExtractData<FEPolarFluidMaterialPoint>();
                    w += ppt->m_gf - ppt->m_hf;
                }
                w /= se->GaussPoints();
                for (int j=0; j<el.Nodes(); ++j)
                    gNodes[el.m_lnode[j]].push_back(h+w);
            }
            else break;
        }
    }
    
    //For each node, average the nodal g
    for (int i=0; i<ps->Nodes(); ++i)
    {
        vec3d g(0,0,0);
        for (int j = 0; j < gNodes[i].size(); ++j)
            g += gNodes[i][j];
        g /= gNodes[i].size();
            
        // store value for now
        m_g[i] = g;
    }
 
    FEPrescribedSurface::Update();
}

//-----------------------------------------------------------------------------
void FEPolarFluidRelativeAngularVelocity::GetNodalValues(int nodelid, std::vector<double>& val)
{
    val[0] = m_g[nodelid].x;
    val[1] = m_g[nodelid].y;
    val[2] = m_g[nodelid].z;
}

//-----------------------------------------------------------------------------
// copy data from another class
void FEPolarFluidRelativeAngularVelocity::CopyFrom(FEBoundaryCondition* pbc)
{
    // TODO: implement this
    assert(false);
}

//-----------------------------------------------------------------------------
//! serialization
void FEPolarFluidRelativeAngularVelocity::Serialize(DumpStream& ar)
{
    FEPrescribedSurface::Serialize(ar);
    ar & m_g;
    if (ar.IsShallow()) return;
    ar & m_dofG;
    ar & m_psurf;
}
