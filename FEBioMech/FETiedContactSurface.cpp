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
#include "FETiedContactSurface.h"
#include "FECore/FEShellDomain.h"
#include "FECore/vector.h"
#include <FECore/FEMesh.h>

FETiedContactSurface::Data::Data()
{
	m_vgap = vec3d(0, 0, 0);
	m_rs = vec2d(0, 0);
	m_Lm = vec3d(0, 0, 0);
	m_Tc = vec3d(0, 0, 0);
	m_off = 0.0;
}

void FETiedContactSurface::Data::Serialize(DumpStream& ar)
{
	FEContactMaterialPoint::Serialize(ar);
	ar & m_vgap;
	ar & m_rs;
	ar & m_Lm;
	ar & m_Tc;
	ar & m_off;
}

//-----------------------------------------------------------------------------
FETiedContactSurface::FETiedContactSurface(FEModel* pfem) : FEContactSurface(pfem)
{ 
	m_boffset = false; 
}

//-----------------------------------------------------------------------------
//! Creates a surface for use with a sliding interface. All surface data
//! structures are allocated.
//! Note that it is assumed that the element array is already created
//! and initialized.

bool FETiedContactSurface::Init()
{
	// always intialize base class first!
	if (FEContactSurface::Init() == false) return false;

	// get the number of nodes
	int nn = Nodes();

	// allocate other surface data
	m_data.resize(nn);

	// we calculate the gap offset values
	// This value is used to take the shell thickness into account
	if (m_boffset)
	{
		FEMesh& m = *m_pMesh;
		vector<double> tag(m.Nodes());
		zero(tag);
		for (int nd = 0; nd<m.Domains(); ++nd)
		{
			FEShellDomain* psd = dynamic_cast<FEShellDomain*>(&m.Domain(nd));
			if (psd)
			{
				for (int i = 0; i<psd->Elements(); ++i)
				{
					FEShellElement& el = psd->Element(i);
					int n = el.Nodes();
					for (int j = 0; j<n; ++j) tag[el.m_node[j]] = 0.5*el.m_h0[j];
				}
			}
		}
		for (int i = 0; i<nn; ++i) m_data[i].m_off = tag[NodeIndex(i)];
	}

	return true;
}

//-----------------------------------------------------------------------------
void FETiedContactSurface::Serialize(DumpStream &ar)
{
	FEContactSurface::Serialize(ar);
	ar & m_data;
}

//-----------------------------------------------------------------------------
void FETiedContactSurface::GetContactTraction(int nface, vec3d& pt)
{
    FESurfaceElement& el = Element(nface);
    int ne = el.Nodes();
    pt = vec3d(0,0,0);
    for (int k=0; k<ne; ++k) pt += m_data[el.m_lnode[k]].m_Tc;
    pt /= ne;
}

//-----------------------------------------------------------------------------
void FETiedContactSurface::GetNodalContactPressure(int nface, double* pn)
{
	FESurfaceElement& f = Element(nface);
	int ne = f.Nodes();
	for (int j= 0; j< ne; ++j) pn[j] = m_data[f.m_lnode[j]].m_Tc.norm();
}

//-----------------------------------------------------------------------------
void FETiedContactSurface::GetNodalContactTraction(int nface, vec3d* tn)
{
	FESurfaceElement& f = Element(nface);
	int ne = f.Nodes();
	for (int j= 0; j< ne; ++j) 
	{
		tn[j] = m_data[f.m_lnode[j]].m_Tc;
	}
}
