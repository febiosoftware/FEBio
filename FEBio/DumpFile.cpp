// DumpFile.cpp: implementation of the DumpFile class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "DumpFile.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

DumpFile::DumpFile()
{
	m_fp = 0;
}

DumpFile::~DumpFile()
{
	Close();
}

bool DumpFile::Open(const char* szfile)
{
	m_fp = fopen(szfile, "rb");
	if (m_fp == 0) return false;

	m_bsave = false;

	return true;
}

bool DumpFile::Create(const char* szfile)
{
	m_fp = fopen(szfile, "wb");
	if (m_fp == 0) return false;

	m_bsave = true;

	return true;
}

bool DumpFile::Append(const char* szfile)
{
	m_fp = fopen(szfile, "a+b");
	if (m_fp == 0) return false;

	m_bsave = true;

	return true;
}

void DumpFile::Close()
{
	if (m_fp) fclose(m_fp); 
	m_fp = 0;
}
