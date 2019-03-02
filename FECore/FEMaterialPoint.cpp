#include "stdafx.h"
#include "FEMaterialPoint.h"
#include <string.h>

FEMaterialPoint::FEMaterialPoint(FEMaterialPoint* ppt)
{
	m_pPrev = 0;
	m_pNext = ppt;
	m_elem = 0;
	if (ppt) ppt->m_pPrev = this;
}

FEMaterialPoint::~FEMaterialPoint()
{ 
	if (m_pNext) delete m_pNext;
	m_pNext = m_pPrev = 0;
}

void FEMaterialPoint::SetPrev(FEMaterialPoint* pt)
{
	m_pPrev = pt;
}

// TODO: What if the next pointer is already assigned?
void FEMaterialPoint::SetNext(FEMaterialPoint* pt)
{
	m_pNext = pt;
	pt->m_pPrev = this;
}

void FEMaterialPoint::Init()
{
	if (m_pNext) m_pNext->Init();
}

void FEMaterialPoint::Update(const FETimeInfo& timeInfo)
{
	if (m_pNext) m_pNext->Update(timeInfo);
}

void FEMaterialPoint::Serialize(DumpStream& ar)
{
	if (m_pNext) m_pNext->Serialize(ar);
}

size_t FEMaterialPoint::memsize()
{
	return sizeof(FEMaterialPoint);
}

//-----------------------------------------------------------------------------
FEMaterialPointArray::FEMaterialPointArray(FEMaterialPoint* ppt) : FEMaterialPoint(ppt)
{
	
}

//-----------------------------------------------------------------------------
void FEMaterialPointArray::AddMaterialPoint(FEMaterialPoint* pt)
{
	m_mp.push_back(pt);
	pt->SetPrev(this);
}


//-----------------------------------------------------------------------------
void FEMaterialPointArray::Init()
{
	for (int i = 0; i<(int)m_mp.size(); ++i) m_mp[i]->Init();
}

//-----------------------------------------------------------------------------
void FEMaterialPointArray::Serialize(DumpStream& ar)
{
	for (int i = 0; i<(int)m_mp.size(); ++i) m_mp[i]->Serialize(ar);
}

//-----------------------------------------------------------------------------
void FEMaterialPointArray::Update(const FETimeInfo& timeInfo)
{
	FEMaterialPoint::Update(timeInfo);
	for (int i = 0; i<(int)m_mp.size(); ++i) m_mp[i]->Update(timeInfo);
}
