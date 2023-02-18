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



#include "stdafx.h"
#include "FEMultiphasicFluidPressureLoad.h"
#include "FESoluteInterface.h"
#include "FEOsmoticCoefficient.h"
#include <FECore/FECoreKernel.h>
#include <FECore/FEModel.h>
#include <FECore/FEMaterial.h>

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEMultiphasicFluidPressureLoad, FESurfaceLoad)
    ADD_PARAMETER(m_p0    , "pressure");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEMultiphasicFluidPressureLoad::FEMultiphasicFluidPressureLoad(FEModel* pfem) : FESurfaceLoad(pfem)
{
    m_p0 = 0;
    
    m_dofP = pfem->GetDOFIndex("p");
    
    m_dof.Clear();
    m_dof.AddDof(m_dofP);
}

//-----------------------------------------------------------------------------
bool FEMultiphasicFluidPressureLoad::Init()
{
    if (FESurfaceLoad::Init() == false) return false;
    m_Rgas = GetFEModel()->GetGlobalConstant("R");
    m_Tabs = GetFEModel()->GetGlobalConstant("T");

    return true;
}

//-----------------------------------------------------------------------------
//! Activate the degrees of freedom for this BC
void FEMultiphasicFluidPressureLoad::Activate()
{
    FESurface* ps = &GetSurface();
    
    for (int i=0; i<ps->Nodes(); ++i)
    {
        FENode& node = ps->Node(i);
        // mark node as having prescribed DOF
        node.set_bc(m_dofP, DOF_PRESCRIBED);
    }
}

//-----------------------------------------------------------------------------
//! Evaluate and prescribe the resistance pressure
void FEMultiphasicFluidPressureLoad::Update()
{
    // prescribe this pressure at the nodes
    FESurface* ps = &GetSurface();
    int N = ps->Nodes();
    vector<double> posm(N,0);
    vector<int> nosm(N,0);
    for (int i=0; i<ps->Elements(); ++i) {
        FESurfaceElement& el = m_psurf->Element(i);
        FEElement* pe = el.m_elem[0];
        FEMaterial* pm = GetFEModel()->GetMaterial(pe->GetMatID());
        FESoluteInterface* psi = pm->ExtractProperty<FESoluteInterface>();
        if (psi) {
            // calculate average osmotic contribution to pressure in this element
            double nint = pe->GaussPoints();
            double p = 0;
            for (int n=0; n<nint; ++n) {
                FEMaterialPoint& mp = *pe->GetMaterialPoint(n);
                double osm = psi->GetOsmolarity(mp);
                double Phi = psi->GetOsmoticCoefficient()->OsmoticCoefficient(mp);
                p -= m_Rgas*m_Tabs*Phi*osm;
            }
            p /= nint;
            for (int i=0; i<el.Nodes(); ++i) {
                int inode = el.m_lnode[i];
                posm[inode] += p;
                ++nosm[inode];
            }
        }
    }

    // now prescribe the effective fluid pressure at surface nodes
    for (int i=0; i<N; ++i)
    {
        if (nosm[i] > 0) posm[i] /= nosm[i];
        if (ps->Node(i).m_ID[m_dofP] < -1)
        {
            FENode& node = ps->Node(i);
            // set nodal value
            node.set(m_dofP, m_p0 + posm[i]);
        }
    }
    
    GetFEModel()->SetMeshUpdateFlag(true);
}

//-----------------------------------------------------------------------------
//! serialization
void FEMultiphasicFluidPressureLoad::Serialize(DumpStream& ar)
{
    FESurfaceLoad::Serialize(ar);
}
