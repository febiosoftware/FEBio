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
#include "FESurfaceAttractionBodyForce.h"
#include "FEElasticMaterial.h"
#include <FECore/FEMesh.h>
#include <FECore/FEClosestPointProjection.h>
#include <FECore/log.h>

BEGIN_FECORE_CLASS(FESurfaceAttractionBodyForce, FEBodyForce);
    ADD_PARAMETER(m_blt     , "blt"          );
    ADD_PARAMETER(m_bsf     , "bsf"          );
    ADD_PARAMETER(m_stol    , "search_tol"   );
    ADD_PARAMETER(m_sradius , "search_radius");

	ADD_PROPERTY(m_s, "surface", FEProperty::Reference);
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FESurfaceAttractionBodyForce::FESurfaceAttractionBodyForce(FEModel* pfem) : FEBodyForce(pfem)
{
    m_blt = 1;
    m_bsf = 0;
    m_stol = 0.01;
    m_sradius = 0;      // no search radius limitation

	m_s = nullptr;
}

//-----------------------------------------------------------------------------
//! initialize
bool FESurfaceAttractionBodyForce::Init()
{
	if (m_s == nullptr) return false;

    if (!m_s->Init()) return false;
    FEBodyLoad::Init();

	FEClosestPointProjection cpp(*m_s);
	cpp.SetTolerance(m_stol);
	cpp.SetSearchRadius(m_sradius);
	cpp.HandleSpecialCases(true);
	cpp.Init();

    // allocate projection point vector array
    int nel = GetMesh().Elements();
    m_q.resize(nel);
    
    const char s[] = "%%";

    feLog("Initializing projections to attractive surface...\n");
    
    // evaluate projection of integration points to attractive surface
    for (int i=0; i<Domains(); ++i) {
        FEDomain* dom = Domain(i);
        int nde = dom->Elements();
        feLog("Domain %d: %d elements\n",i+1,nde);
        int pcp = 0;
        feLog("\r%d %s done",pcp,s);
        for (int j=0; j<nde; ++j) {
            int pc = 100*(j+1)/nde;
            if (pc > pcp) {
                feLog("\r%d %s done",pc,s);
                pcp = pc;
            }
            FEElement& el = dom->ElementRef(j);
            int nint = el.GaussPoints();
            int eid = el.GetID() - 1;
            m_q[eid].resize(nint);
            for (int k=0; k<nint; ++k) {
                FEMaterialPoint* mp = el.GetMaterialPoint(k);
                vec3d x = mp->m_r0;

				vec2d rs(0, 0);
				FESurfaceElement* pme = nullptr;
				pme = cpp.Project(x, m_q[eid][k], rs);

                if (pme == nullptr)
                    m_q[eid][k] = x + vec3d(100*m_blt,100*m_blt,100*m_blt);
            }
        }
    }
    
    return true;
}

//-----------------------------------------------------------------------------
vec3d FESurfaceAttractionBodyForce::force(FEMaterialPoint& mp)
{
    // get element number for this material point
    int eid = mp.m_elem->GetID() - 1;
    
    vec3d q = m_q[eid][mp.m_index];
    
    // initialize net force
    vec3d f(0,0,0);

    // calculate net force
    vec3d g = mp.m_r0 - q;
    double r = g.unit();
    f = g*(m_bsf*exp(-r/m_blt));

    return f;
}

mat3ds FESurfaceAttractionBodyForce::stiffness(FEMaterialPoint& mp)
{
    return mat3ds(0,0,0,0,0,0);
}
