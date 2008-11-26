// Archive.cpp: implementation of the Archive class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "Archive.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Archive::Archive()
{
	m_fp = 0;
}

Archive::~Archive()
{
	Close();
}

bool Archive::Open(const char* szfile)
{
	m_fp = fopen(szfile, "rb");
	if (m_fp == 0) return false;

	m_bsave = false;

	return true;
}

bool Archive::Create(const char* szfile)
{
	m_fp = fopen(szfile, "wb");
	if (m_fp == 0) return false;

	m_bsave = true;

	return true;
}

bool Archive::Append(const char* szfile)
{
	m_fp = fopen(szfile, "a+b");
	if (m_fp == 0) return false;

	m_bsave = true;

	return true;
}

void Archive::Close()
{
	if (m_fp) fclose(m_fp); 
	m_fp = 0;
}
