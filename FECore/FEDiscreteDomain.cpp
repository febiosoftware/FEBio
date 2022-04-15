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
#include "FEDiscreteDomain.h"
#include "FEMesh.h"
#include "FEMaterial.h"

//-----------------------------------------------------------------------------
bool FEDiscreteDomain::Create(int nelems, FE_Element_Spec espec)
{ 
	m_Elem.resize(nelems); 
	for (int i = 0; i < nelems; ++i)
	{
		FEDiscreteElement& el = m_Elem[i];
		el.SetLocalID(i);
		el.SetMeshPartition(this);
	}

	if (espec.etype != FE_ELEM_INVALID_TYPE)
		for (int i=0; i<nelems; ++i) m_Elem[i].SetType(espec.etype);

	return true;
}

//-----------------------------------------------------------------------------
void FEDiscreteDomain::Reset()
{
	for (auto& el : m_Elem) el.setActive();
}

//-----------------------------------------------------------------------------
bool FEDiscreteDomain::Init()
{
	if (FEDomain::Init() == false) return false;

	FEMaterial* pmat = GetMaterial();
	if (pmat) SetMatID(pmat->GetID());

	return true;
}

//-----------------------------------------------------------------------------
void FEDiscreteDomain::CopyFrom(FEMeshPartition* pd)
{
	FEDomain::CopyFrom(pd);
	FEDiscreteDomain* psd = dynamic_cast<FEDiscreteDomain*>(pd);
	m_Elem = psd->m_Elem;
	ForEachElement([=](FEElement& el) { el.SetMeshPartition(this); });
}

//-----------------------------------------------------------------------------
void FEDiscreteDomain::AddElement(int eid, int n[2])
{
	FEDiscreteElement el;
	el.SetType(FE_DISCRETE);
	el.m_node[0] = n[0];
	el.m_node[1] = n[1];
	el.SetID(eid);
	m_Elem.push_back(el);
}
