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
#include "FEMatchingOsmoticCoefficientBC.h"
#include "FEBioMix.h"
#include <FECore/FEModel.h>

//=============================================================================
BEGIN_FECORE_CLASS(FEMatchingOsmoticCoefficientBC, FEPrescribedSurface)
    ADD_PARAMETER(m_ambp , "ambient_pressure");
    ADD_PARAMETER(m_ambc , "ambient_osmolarity");
    ADD_PARAMETER(m_bshellb , "shell_bottom");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FEMatchingOsmoticCoefficientBC::FEMatchingOsmoticCoefficientBC(FEModel* pfem) : FEPrescribedSurface(pfem)
{
    m_ambp = m_ambc = 0.0;
    m_bshellb = false;
    
    m_dofP = pfem->GetDOFIndex("p");
    m_dofQ = pfem->GetDOFIndex("q");
}

//-----------------------------------------------------------------------------
//! Activate the degrees of freedom for this BC
void FEMatchingOsmoticCoefficientBC::Activate()
{
    FESurface* ps = GetSurface();
    
    FEModel* fem = GetFEModel();

    for (int i=0; i<ps->Elements(); ++i) {
        FESurfaceElement& el = ps->Element(i);
        // get the element connected to this surface
        FEElement* elem = el.m_elem[0];
        FEMaterial* pm = fem->GetMaterial(elem->GetMatID());
        // get the multihasic material for this element
        FEMultiphasic* pmp = dynamic_cast<FEMultiphasic*>(pm);
        // if the material is multiphasic, set prescribable fluid dofs
        if (pmp) {
            for (int j=0; j<el.Nodes(); ++j) {
                FENode& node = ps->Node(el.m_lnode[j]);
                // mark node as having prescribed DOF
                if (!m_bshellb) node.set_bc(m_dofP, DOF_PRESCRIBED);
                else node.set_bc(m_dofQ, DOF_PRESCRIBED);
            }
        }
    }

    // NOTE: Hmmm, not sure about this one. Maybe skip immediate base class?
//    FEPrescribedSurface::Activate();
    FESurfaceBC::Activate();
}

//-----------------------------------------------------------------------------
//! Evaluate and prescribe the ambient effective fluid pressure, using the multiphasic osmotic coefficient
void FEMatchingOsmoticCoefficientBC::Update()
{
    FESurface* ps = GetSurface();
    
    FEModel* fem = GetFEModel();
    
    vector<double> pn(ps->Nodes(),0);
    vector<int> cnt(ps->Nodes(),0);

    for (int i=0; i<ps->Elements(); ++i) {
        FESurfaceElement& el = ps->Element(i);
        // get the element connected to this surface
        FEElement* elem = el.m_elem[0];
        FEMaterial* pm = fem->GetMaterial(elem->GetMatID());
        // get the multihasic material for this element
        FEMultiphasic* pmp = dynamic_cast<FEMultiphasic*>(pm);
        double Phi = 0;
        if (pmp) {
            // get average osmotic coefficient in this multiphasic element
            int nint = elem->GaussPoints();
            for (int n=0; n<nint; ++n) {
                FEMaterialPoint& mp = *elem->GetMaterialPoint(n);
                Phi += pmp->GetOsmoticCoefficient()->OsmoticCoefficient(mp);
            }
            Phi /= nint;
            // evaluate effective fluid pressure on the boundary
            double pe = m_ambp - pmp->m_Rgas*pmp->m_Tabs*Phi*m_ambc;
            // project to nodes of that face
            for (int j=0; j<el.Nodes(); ++j) {
                int n = el.m_lnode[j];
                pn[n] += pe;
                ++cnt[n];
            }
        }
    }

    // prescribe the projected nodal values of the effective fluid pressure
    for (int i=0; i<ps->Nodes(); ++i)
    {
        if (ps->Node(i).m_ID[m_dofP] < -1)
        {
            assert(cnt[i]);
            FENode& node = ps->Node(i);
            // prescribe effective pressure at this node
            double pe = pn[i]/cnt[i];
            if (!m_bshellb) node.set(m_dofP, pe);
            else node.set(m_dofQ, pe);
        }
    }
}

//-----------------------------------------------------------------------------
void FEMatchingOsmoticCoefficientBC::GetNodalValues(int nodelid, std::vector<double>& val)
{
    // TODO: implement this
    assert(false);
}

//-----------------------------------------------------------------------------
void FEMatchingOsmoticCoefficientBC::CopyFrom(FEBoundaryCondition* pbc)
{
    // TODO: implement this
    assert(false);
}

//-----------------------------------------------------------------------------
//! serialization
void FEMatchingOsmoticCoefficientBC::Serialize(DumpStream& ar)
{
    FEPrescribedSurface::Serialize(ar);
    if (ar.IsShallow() == false)
    {
        ar& m_dofP;
        ar& m_dofQ;
    }
}
