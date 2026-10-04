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
    ar & m_dg;
    ar & m_Lmd;
    ar & m_Lmt;
    ar & m_epsn;
    ar & m_epsp;
    ar & m_p1;
    ar & m_nu;
    ar & m_s1;
    ar & m_tr;
    ar & m_rs;
    ar & m_rsp;
    ar & m_bstick;
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
void FEBiphasicContactSurface::GetGPLocalFLS(int nface, double* fls, double pamb)
{
    FESurfaceElement& el = Element(nface);

    const int nint = el.GaussPoints();
    const int neln = el.Nodes();

    // NOTE: the output array is indexed by integration point, so it must be
    //       zeroed over nint entries.  This used to zero el.Nodes() entries,
    //       which left fls[neln..nint-1] uninitialized on any facet with more
    //       integration points than nodes (QUAD8G9, for example).
    for (int i=0; i<nint; ++i) fls[i] = 0.0;

    FEElement* e = el.m_elem[0].pe;
    FESolidElement* se = dynamic_cast<FESolidElement*>(e);
    if (se == nullptr) return;

    // average the effective stress and the fluid pressure over the parent
    // solid element
    mat3ds s; s.zero();
    double p = 0;
    for (int i=0; i<se->GaussPoints(); ++i) {
        FEMaterialPoint* mp = se->GetMaterialPoint(i);
        FEElasticMaterialPoint* ep = mp->ExtractData<FEElasticMaterialPoint>();
        FEBiphasicMaterialPoint* bp = mp->ExtractData<FEBiphasicMaterialPoint>();
        if (ep) s += ep->m_s;
        if (bp) p += bp->m_p;
    }
    s /= se->GaussPoints();
    p /= se->GaussPoints();

    // account for ambient pressure
    p -= pamb;

    // Evaluate a normal at each node of the facet, averaged from the facet's
    // integration points and weighted by that node's shape function, so the
    // integration points nearest the node dominate.
    //
    // NOTE: FESurface::SurfaceNormal(el, n) takes an INTEGRATION POINT index,
    //       not a node index.  This function used to call it as
    //       SurfaceNormal(el, j) with j a node index, which indexes the shape
    //       function derivative tables (el.Gr(n), el.Gs(n)) out of range on
    //       every facet type with fewer integration points than nodes --
    //       including TRI3, which FESurface::Create() builds as FE_TRI3G1:
    //       one integration point, three nodes.  That is an out-of-bounds read
    //       followed by a dereference of whatever it returns.
    vec3d nn[FEElement::MAX_NODES];
    for (int j=0; j<neln; ++j) nn[j] = vec3d(0,0,0);
    for (int i=0; i<nint; ++i) {
        vec3d ni = SurfaceNormal(el, i);
        double* H = el.H(i);
        for (int j=0; j<neln; ++j) nn[j] += ni*H[j];
    }

    // evaluate the FLS at the nodes
    double flsn[FEElement::MAX_NODES];
    for (int j=0; j<neln; ++j) {
        // vec3d::unit() leaves a zero vector alone, which would give tn = 0
        // and hence fls = 0 below -- the same result as a degenerate facet.
        nn[j].unit();
        double tn = nn[j]*(s*nn[j]);
        flsn[j] = (tn != 0) ? -p/tn : 0;
    }

    // interpolate to the integration points of the facet
    for (int i=0; i<nint; ++i) {
        double* H = el.H(i);
        for (int j=0; j<neln; ++j) fls[i] += flsn[j]*H[j];
    }
}
