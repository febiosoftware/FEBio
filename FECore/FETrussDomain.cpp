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
#include "FETrussDomain.h"
#include "FEMesh.h"

//-----------------------------------------------------------------------------
FETrussDomain::FETrussDomain(FEModel* fem) : FEBeamDomain(fem)
{
}

//-----------------------------------------------------------------------------
bool FETrussDomain::Create(int nsize, FE_Element_Spec espec)
{
	m_Elem.resize(nsize);
	for (int i = 0; i < nsize; ++i)
	{
		FETrussElement& el = m_Elem[i];
		el.SetLocalID(i);
		el.SetMeshPartition(this);
	}

	if (espec.etype != FE_ELEM_INVALID_TYPE)
		for (int i=0; i<nsize; ++i) m_Elem[i].SetType(espec.etype);

	return true;
}

//-----------------------------------------------------------------------------
vec3d FETrussDomain::TrussNormal(FETrussElement& el)
{
	vec3d rt[2];
	rt[0] = m_pMesh->Node(el.m_node[0]).m_rt;
	rt[1] = m_pMesh->Node(el.m_node[1]).m_rt;

	vec3d a = rt[0];
	vec3d b = rt[1];
	vec3d n = b - a;
	n.unit();
	return n;
}

//-----------------------------------------------------------------------------
void FETrussDomain::ForEachTrussElement(std::function<void(FETrussElement& el)> f)
{
	int N = Elements();
	for (int i = 0; i < N; ++i)
	{
		FETrussElement& el = m_Elem[i];
		f(el);
	}
}
