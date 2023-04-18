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
#include "FEBiphasicContactSurface.h"
#include "FEBiphasic.h"
#include <FECore/FEModel.h>

void FEBiphasicContactPoint::Serialize(DumpStream& ar)
{
    FEContactMaterialPoint::Serialize(ar);
    ar & m_Lmp & m_pg & m_mueff & m_fls;
}

//-----------------------------------------------------------------------------
FEBiphasicContactSurface::FEBiphasicContactSurface(FEModel* pfem) : FEContactSurface(pfem)
{
	m_dofP = -1;
}

//-----------------------------------------------------------------------------
FEBiphasicContactSurface::~FEBiphasicContactSurface()
{
}

//-----------------------------------------------------------------------------
bool FEBiphasicContactSurface::Init()
{
	// I want to use the FEModel class for this, but don't know how
	DOFS& dofs = GetFEModel()->GetDOFS();
	m_dofP = dofs.GetDOF("p");
	return FEContactSurface::Init();
}

//-----------------------------------------------------------------------------
//! serialization
void FEBiphasicContactSurface::Serialize(DumpStream& ar)
{
	FEContactSurface::Serialize(ar);
	if (ar.IsShallow() == false) ar & m_dofP;
}

//-----------------------------------------------------------------------------
vec3d FEBiphasicContactSurface::GetFluidForce()
{
	assert(false);
    return vec3d(0,0,0);
}

//-----------------------------------------------------------------------------
double FEBiphasicContactSurface::GetFluidLoadSupport()
{
    int n, i;
    
    // initialize contact force
    double FLS = 0;
    double A = 0;
    
    // loop over all elements of the surface
    for (n=0; n<Elements(); ++n)
    {
        FESurfaceElement& el = Element(n);
        // evaluate the fluid force for that element
        for (i=0; i<el.GaussPoints(); ++i)
        {
            FEBiphasicContactPoint *cp = dynamic_cast<FEBiphasicContactPoint*>(el.GetMaterialPoint(i));
            if (cp) {
                double w = el.GaussWeights()[i];
                // get the base vectors
                vec3d g[2];
                CoBaseVectors(el, i, g);
                // normal (magnitude = area)
                vec3d n = g[0] ^ g[1];
                double da = n.norm();
                FLS += cp->m_fls*w*da;
            }
        }
    }
    
    A = GetContactArea();
    
    return (A > 0) ? FLS/A : 0;
}

//-----------------------------------------------------------------------------
void FEBiphasicContactSurface::GetMuEffective(int nface, double& pg)
{
    pg = 0;
}

//-----------------------------------------------------------------------------
void FEBiphasicContactSurface::GetLocalFLS(int nface, double& pg)
{
    pg = 0;
}

//-----------------------------------------------------------------------------
void FEBiphasicContactSurface::UnpackLM(FEElement& el, vector<int>& lm)
{
	int N = el.Nodes();
	lm.assign(N*4, -1);

	// pack the equation numbers
	for (int i=0; i<N; ++i)
	{
		int n = el.m_node[i];

		FENode& node = m_pMesh->Node(n);
		vector<int>& id = node.m_ID;

		// first the displacement dofs
		lm[3*i  ] = id[m_dofX];
		lm[3*i+1] = id[m_dofY];
		lm[3*i+2] = id[m_dofZ];

		// now the pressure dofs
		if (m_dofP >= 0) lm[3*N+i] = id[m_dofP];
	}
}
//-----------------------------------------------------------------------------
// Evaluate the local fluid load support projected from the element to the surface Gauss points
void FEBiphasicContactSurface::GetGPLocalFLS(int nface, double* pt, double pamb)
{
    FESurfaceElement& el = Element(nface);
    FEElement* e = el.m_elem[0];
    FESolidElement* se = dynamic_cast<FESolidElement*>(e);
    if (se) {
        mat3ds s; s.zero();
        double p = 0;
        for (int i=0; i<se->GaussPoints(); ++i) {
            FEMaterialPoint* pt = se->GetMaterialPoint(i);
            FEElasticMaterialPoint* ep = pt->ExtractData<FEElasticMaterialPoint>();
            FEBiphasicMaterialPoint* bp = pt->ExtractData<FEBiphasicMaterialPoint>();
            if (ep) s += ep->m_s;
            if (bp) p += bp->m_p;
        }
        s /= se->GaussPoints();
        p /= se->GaussPoints();
        // account for ambient pressure
        p -= pamb;
        // evaluate FLS at integration points of that face
        for (int i=0; i<el.GaussPoints(); ++i) {
            double *H = el.H(i);
            pt[i] = 0;
            for (int j=0; j<el.Nodes(); ++j) {
                vec3d n = SurfaceNormal(el, j);
                double tn = n*(s*n);
                double fls = (tn != 0) ? -p/tn : 0;
                pt[i] += fls*H[j];
            }
        }
    }
    else
        for (int i=0; i<el.Nodes(); ++i) pt[i] = 0;
}
