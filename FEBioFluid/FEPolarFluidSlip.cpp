/*This file is part of the FEBio source code and is licensed under the MIT license
 listed below.
 
 See Copyright-FEBio.txt for details.
 
 Copyright (c) 2020 University of Utah, The Trustees of Columbia University in
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
#include "FEPolarFluidSlip.h"
#include "FEBioPolarFluid.h"
#include <FECore/FESurface.h>
#include <FECore/FEModel.h>
#include "FEPolarFluid.h"

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEPolarFluidSlip, FEPrescribedSurface)
    ADD_PARAMETER(m_m0, "m0")->setUnits(UNIT_STIFFNESS)->setLongName("couple traction threshold");
    ADD_PARAMETER(m_ksi, "ksi")->setUnits("t/M")->setLongName("slope");
    ADD_PARAMETER(m_brel, "relative");
END_FECORE_CLASS();


//-----------------------------------------------------------------------------
//! constructor
FEPolarFluidSlip::FEPolarFluidSlip(FEModel* fem) : FEPrescribedSurface(fem), m_dofW(fem), m_dofG(fem)
{
    m_m0 = 0.0;
    m_ksi = 0.0;
    m_psurf = nullptr;

    m_dofW.AddVariable(FEBioPolarFluid::GetVariableName(FEBioPolarFluid::RELATIVE_FLUID_VELOCITY));
    m_dofG.AddVariable(FEBioPolarFluid::GetVariableName(FEBioPolarFluid::FLUID_ANGULAR_VELOCITY));
    m_dof.AddDofs(m_dofW);
    m_dof.AddDofs(m_dofG);
}

//-----------------------------------------------------------------------------
//! initialize
bool FEPolarFluidSlip::Init()
{
    
    SetRelativeFlag(m_brel);
    if (FEPrescribedSurface::Init() == false) return false;
    
    m_psurf = GetSurface();
    
    m_v.assign(GetSurface()->Nodes(), vec3d(0,0,0));
    m_g.assign(GetSurface()->Nodes(), vec3d(0,0,0));
    
    // do an initial Update so that the velocities and angular velocities are set properly
    // at the very first time step
    Update();
    
    return true;
}

//-----------------------------------------------------------------------------
void FEPolarFluidSlip::UpdateAngularVelocity()
{
    int N = m_psurf->Nodes();
    std::vector<vector<vec3d>> vNodes(N, vector<vec3d>());
    std::vector<vector<vec3d>> gNodes(N, vector<vec3d>());
    
    //Project g value from int points to nodes on surface
    for (int i = 0; i < m_psurf->Elements(); ++i)
    {
        FESurfaceElement& el = m_psurf->Element(i);
        FEElement* pe = el.m_elem[0];
        vec3d nu = m_psurf->SurfaceNormal(el,0,0);
        if (pe) {
            // get the material
            FEMaterial* pm = GetFEModel()->GetMaterial(pe->GetMatID());
            FEFluidMaterial* pfluid = pm->ExtractProperty<FEFluidMaterial>();
            
            if (!pfluid) {
                pe = el.m_elem[1];
                if (pe) pfluid = GetFEModel()->GetMaterial(pe->GetMatID())->ExtractProperty<FEFluidMaterial>();
            }
            
            // see if this is a fluid element
            if (pfluid) {
                FEPolarFluidMaterial* polar = pfluid->ExtractProperty<FEPolarFluidMaterial>();
                if (polar) {
                    FEPolarFluid* pf = polar->ExtractProperty<FEPolarFluid>();
                    // evaluate the average stress in this element
                    int nint = pe->GaussPoints();
                    mat3d M(mat3dd(0));
                    vec3d w(0,0,0);
                    for (int n=0; n<nint; ++n)
                    {
                        FEMaterialPoint& mp = *pe->GetMaterialPoint(n);
                        FEFluidMaterialPoint& pf = *(mp.ExtractData<FEFluidMaterialPoint>());
                        mat3da W(pf.m_Lf);
                        w += W.vec();
                        FEPolarFluidMaterialPoint& pt = *(mp.ExtractData<FEPolarFluidMaterialPoint>());
                        M += pt.m_M;
                    }
                    M /= nint;
                    w /= nint;
                    // value of couple traction on this surface element
                    vec3d mt = M*nu;
                    vec3d v(0,0,0);
                    vec3d g(0,0,0);
                    // if above threshold, calculate tangential angular velocity
                    double mtlen = mt.unit();
                    if (mtlen > m_m0) {
                        g = mt*((mtlen - m_m0)*m_ksi);
                        v = g ^ (nu*pf->m_kg); // evaluate slip velocity from pure rolling
                    }
                    for (int j=0; j<el.Nodes(); ++j) {
                        gNodes[el.m_lnode[j]].push_back(g);
                        vNodes[el.m_lnode[j]].push_back(v);
                    }
                }
            }
        }
    }
    
    //For each node, average the nodal g
    for (int i = 0; i < m_psurf->Nodes(); ++i)
    {
        vec3d v(0,0,0);
        vec3d g(0,0,0);
        for (int j = 0; j < gNodes[i].size(); ++j) {
            v += vNodes[i][j];
            g += gNodes[i][j];
        }
        if (vNodes[i].size()) v /= vNodes[i].size();
        if (gNodes[i].size()) g /= gNodes[i].size();

        // store value for now
        m_v[i] = v;
        m_g[i] = g;
    }
}

//-----------------------------------------------------------------------------
//! Evaluate and prescribe the resistance pressure
void FEPolarFluidSlip::UpdateModel() { Update(); }
void FEPolarFluidSlip::Update()
{
    UpdateAngularVelocity();
    FEPrescribedSurface::Update();
    
    // TODO: Is this necessary?
    GetFEModel()->SetMeshUpdateFlag(true);
}

//-----------------------------------------------------------------------------
void FEPolarFluidSlip::PrepStep(std::vector<double>& ui, bool brel)
{
    UpdateAngularVelocity();
    FEPrescribedSurface::PrepStep(ui, brel);
}

//-----------------------------------------------------------------------------
void FEPolarFluidSlip::GetNodalValues(int nodelid, std::vector<double>& val)
{
    val[0] = m_v[nodelid].x;
    val[1] = m_v[nodelid].y;
    val[2] = m_v[nodelid].z;
    val[3] = m_g[nodelid].x;
    val[4] = m_g[nodelid].y;
    val[5] = m_g[nodelid].z;
}

void FEPolarFluidSlip::CopyFrom(FEBoundaryCondition* pbc)
{
    // TODO: implement
    assert(false);
}

//-----------------------------------------------------------------------------
//! serialization
void FEPolarFluidSlip::Serialize(DumpStream& ar)
{
    FEPrescribedSurface::Serialize(ar);
    ar & m_v;
    ar & m_g;
    if (ar.IsShallow()) return;
    ar & m_dofW;
    ar & m_dofG;
    ar & m_psurf;
}
