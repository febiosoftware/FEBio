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
#include "FEPressureStabilization.h"
#include "FEBiphasic.h"
#include "FECore/FEModel.h"

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEPressureStabilization, FESurfaceLoad)
	ADD_PARAMETER(m_bstab, "stabilize"   );
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FEPressureStabilization::FEPressureStabilization(FEModel* pfem) : FESurfaceLoad(pfem)
{
    m_bstab = true;
}

//-----------------------------------------------------------------------------
//! allocate storage
void FEPressureStabilization::SetSurface(FESurface* ps)
{
    FESurfaceLoad::SetSurface(ps);
}

//-----------------------------------------------------------------------------
bool FEPressureStabilization::Init()
{
    FESurface& ps = GetSurface();
    ps.Init();
    
    return true;
}

//-----------------------------------------------------------------------------
void FEPressureStabilization::Activate()
{
    FESurface& ps = GetSurface();
    
    // get the mesh
    FEMesh& m = GetFEModel()->GetMesh();
    
    // loop over all surface elements
    for (int i=0; i<ps.Elements(); ++i)
    {
        // get the surface element
        FESurfaceElement& el = ps.Element(i);
        
        // find the element this face belongs to
        FEElement* pe = el.m_elem[0];
        assert(pe);
        
        // calculate time constant
        double tau = TimeConstant(el, ps);
    }
    
    FESurfaceLoad::Activate();
}

//-----------------------------------------------------------------------------
double FEPressureStabilization::TimeConstant(FESurfaceElement& el, FESurface& s)
{
    // get the mesh
    FEMesh& m = GetFEModel()->GetMesh();
    
    double tau = 0;
    
    // get the element this surface element belongs to
    FEElement* pe = el.m_elem[0];
    if (pe)
    {
        // evaluate element surface normal at parametric center
        vec3d t[2];
        s.CoBaseVectors0(el, 0, 0, t);
        vec3d n = t[0] ^ t[1];
        n.unit();
        
        // get the material
        FEMaterial* pm = GetFEModel()->GetMaterial(pe->GetMatID());
        
        // see if this is a poro-elastic element
        FEBiphasic* biph = dynamic_cast<FEBiphasic*> (pm);
        if (biph)
        {
            // get the area of the surface element
            double A = s.FaceArea(el);
            
            // get the volume of the volume element
            double V = m.ElementVolume(*pe);
            
            // calculate the element thickness
            double h = V/A;
            
            // get a material point
            FEMaterialPoint& mp = *pe->GetMaterialPoint(0);
            FEElasticMaterialPoint& ept = *(mp.ExtractData<FEElasticMaterialPoint>());
            
            // setup the material point
            ept.m_F = mat3dd(1.0);
            ept.m_J = 1;
            ept.m_s.zero();
            
            // get the tangent (stiffness) and it inverse (compliance) at this point
            tens4ds C = biph->Tangent(mp);
            double Ha = n*(vdotTdotv(n, C, n)*n);
            
            // if this is a poroelastic element, then get the permeability tensor
            FEBiphasicMaterialPoint& pt = *(mp.ExtractData<FEBiphasicMaterialPoint>());
            pt.m_p = 0;
            pt.m_w = vec3d(0,0,0);
            
            mat3ds K = biph->Permeability(mp);

            double k = n*(K*n);
            
            tau = h*h/(4*Ha*k);
            
            // if time constant not yet set
            if (biph->m_tau == 0)
                biph->m_tau = tau;  // set it to calculated value
            else
                // pick smallest value
                biph->m_tau = max(biph->m_tau, tau);
        }
    }
    
    return tau;
}

