/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, Columbia University, and others.

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
#include "FEFluidNormalTraction.h"
#include "FEBioFluid.h"
#include <FECore/FESurface.h>
#include <FECore/FEFacetSet.h>
#include <FECore/FEMesh.h>

//=============================================================================
BEGIN_FECORE_CLASS(FEFluidNormalTraction, FESurfaceLoad)
	ADD_PARAMETER(m_traction, "traction");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FEFluidNormalTraction::FEFluidNormalTraction(FEModel* pfem) : FESurfaceLoad(pfem), m_dofW(pfem)
{
	m_traction = 1.0;
    
	m_dofW.AddVariable(FEBioFluid::GetVariableName(FEBioFluid::RELATIVE_FLUID_VELOCITY));
}

//-----------------------------------------------------------------------------
//! allocate storage
void FEFluidNormalTraction::SetSurface(FESurface* ps)
{
    FESurfaceLoad::SetSurface(ps);
	m_traction.SetItemList(ps->GetFacetSet());
}

//-----------------------------------------------------------------------------
void FEFluidNormalTraction::UnpackLM(FEElement& el, vector<int>& lm)
{
    FEMesh& mesh = *GetSurface().GetMesh();
    int N = el.Nodes();
    lm.resize(N*3);
    for (int i=0; i<N; ++i)
    {
        int n = el.m_node[i];
        FENode& node = mesh.Node(n);
        vector<int>& id = node.m_ID;
        
        lm[3*i  ] = id[m_dofW[0]];
        lm[3*i+1] = id[m_dofW[1]];
        lm[3*i+2] = id[m_dofW[2]];
    }
}

//-----------------------------------------------------------------------------
//! Calculate the residual for the traction load
void FEFluidNormalTraction::Residual(FEGlobalVector& R, const FETimeInfo& tp)
{
    vector<double> fe;
    vector<int> elm;
    
    vec3d rt[FEElement::MAX_NODES];
    
    int i, n;
    int N = m_psurf->Elements();
    for (int iel=0; iel<N; ++iel)
    {
        FESurfaceElement& el = m_psurf->Element(iel);
        
        int ndof = 3*el.Nodes();
        fe.resize(ndof);
        
        // nr integration points
        int nint = el.GaussPoints();
        
        // nr of element nodes
        int neln = el.Nodes();
        
        // nodal coordinates
        for (i=0; i<neln; ++i) {
            FENode& node = m_psurf->GetMesh()->Node(el.m_node[i]);
            rt[i] = node.m_rt*tp.alphaf + node.m_rp*(1-tp.alphaf);
        }
        
        double* Gr, *Gs;
        double* N;
        double* w  = el.GaussWeights();
        
        vec3d dxr, dxs;
        
        // repeat over integration points
        zero(fe);
        for (n=0; n<nint; ++n)
        {
			FEMaterialPoint& mp = *el.GetMaterialPoint(n);

            N  = el.H(n);
            Gr = el.Gr(n);
            Gs = el.Gs(n);
            
            // calculate the tangent vectors
            dxr = dxs = vec3d(0,0,0);
            for (i=0; i<neln; ++i)
            {
                dxr.x += Gr[i]*rt[i].x;
                dxr.y += Gr[i]*rt[i].y;
                dxr.z += Gr[i]*rt[i].z;
                
                dxs.x += Gs[i]*rt[i].x;
                dxs.y += Gs[i]*rt[i].y;
                dxs.z += Gs[i]*rt[i].z;
            }

			double tn = m_traction(mp);
            
            vec3d normal = dxr ^ dxs;
            
            for (i=0; i<neln; ++i)
            {
				vec3d f = normal*(tn*w[n]);

                fe[3*i  ] += N[i]*f.x;
                fe[3*i+1] += N[i]*f.y;
                fe[3*i+2] += N[i]*f.z;
            }
        }
        
        // get the element's LM vector and adjust it
        UnpackLM(el, elm);
        
        // add element force vector to global force vector
        R.Assemble(el.m_node, elm, fe);
    }
}
