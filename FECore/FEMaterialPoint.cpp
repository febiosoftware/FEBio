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
#include "FEMaterialPoint.h"
#include "DumpStream.h"
#include <string.h>

FEMaterialPointData::FEMaterialPointData(FEMaterialPointData* ppt)
{
	m_pPrev = 0;
	m_pNext = ppt;
	if (ppt) ppt->m_pPrev = this;
}

FEMaterialPointData::~FEMaterialPointData()
{ 
	if (m_pNext) delete m_pNext;
	m_pNext = m_pPrev = 0;
}

void FEMaterialPointData::SetPrev(FEMaterialPointData* pt)
{
	m_pPrev = pt;
}

void FEMaterialPointData::SetNext(FEMaterialPointData* pt)
{
	m_pNext = pt;
	if (pt) pt->m_pPrev = this;
}

void FEMaterialPointData::Append(FEMaterialPointData* pt)
{
	if (pt == nullptr) return;
	if (m_pNext) m_pNext->Append(pt);
	else SetNext(pt);
}

void FEMaterialPointData::Init()
{
	if (m_pNext) m_pNext->Init();
}

void FEMaterialPointData::Update(const FETimeInfo& timeInfo)
{
	if (m_pNext) m_pNext->Update(timeInfo);
}

void FEMaterialPointData::Serialize(DumpStream& ar)
{
	if (m_pNext) m_pNext->Serialize(ar);
}

//=================================================================================================
FEMaterialPoint::FEMaterialPoint(FEMaterialPointData* data)
{
	m_data = data;
	m_elem = nullptr;
}

FEMaterialPoint::~FEMaterialPoint()
{

}

void FEMaterialPoint::Init()
{
	if (m_data) m_data->Init();
}

FEMaterialPoint* FEMaterialPoint::Copy()
{
	FEMaterialPoint* mp = new FEMaterialPoint(*this);
	if (m_data) mp->m_data = m_data->Copy();
	return mp;
}

void FEMaterialPoint::Update(const FETimeInfo& timeInfo)
{
	if (m_data) m_data->Update(timeInfo);
}

void FEMaterialPoint::Serialize(DumpStream& ar)
{
	if (ar.IsShallow() == false)
	{
		ar & m_r0 & m_J0 & m_Jt;
	}
	if (m_data) m_data->Serialize(ar);
}

void FEMaterialPoint::Append(FEMaterialPointData* pt)
{
	if (pt == nullptr) return;
	assert(m_data);
	if (m_data) m_data->Append(pt);
}

//=================================================================================================

//-----------------------------------------------------------------------------
FEMaterialPointArray::FEMaterialPointArray(FEMaterialPointData* ppt) : FEMaterialPointData(ppt)
{
	
}

//-----------------------------------------------------------------------------
void FEMaterialPointArray::AddMaterialPoint(FEMaterialPoint* pt)
{
	m_mp.push_back(pt);
}

//-----------------------------------------------------------------------------
void FEMaterialPointArray::Init()
{
	for (int i = 0; i<(int)m_mp.size(); ++i) m_mp[i]->Init();
	FEMaterialPointData::Init();
}

//-----------------------------------------------------------------------------
void FEMaterialPointArray::Serialize(DumpStream& ar)
{
	FEMaterialPointData::Serialize(ar);
	for (int i = 0; i<(int)m_mp.size(); ++i) m_mp[i]->Serialize(ar);
}

//-----------------------------------------------------------------------------
void FEMaterialPointArray::Update(const FETimeInfo& timeInfo)
{
	FEMaterialPointData::Update(timeInfo);
	for (int i = 0; i<(int)m_mp.size(); ++i) m_mp[i]->Update(timeInfo);
}
