// DumpFile.cpp: implementation of the DumpFile class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "DumpFile.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

DumpFile::DumpFile(FEModel* pfem)
{
	m_fp = 0;
	m_pfem = pfem;
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

DumpFile& DumpFile::operator << (const char* sz) 
{ 
	int n = (sz ? strlen(sz) : 0); 
	fwrite(&n, sizeof(int), 1, m_fp);
	if (sz) fwrite(sz, sizeof(char), n, m_fp);
	return (*this);
}

DumpFile& DumpFile::operator << (char* sz) 
{ 
	int n = (sz ? strlen(sz) : 0); 
	fwrite(&n, sizeof(int), 1, m_fp);
	if (sz) fwrite(sz, sizeof(char), n, m_fp);
	return (*this);
}

DumpFile& DumpFile::operator << (const double a[3][3])
{
	fwrite(a, sizeof(double), 9, m_fp);
	return (*this);
}

DumpFile& DumpFile::operator << (std::vector<bool>& v)
{
	int n = v.size();
	fwrite(&n, sizeof(int), 1, m_fp);
	for (int i=0; i<n; ++i) { int a = (v[i]?1:0); fwrite(&a, sizeof(int), 1, m_fp); }
	return (*this);
}

DumpFile& DumpFile::operator >> (char* sz) 
{ 
	int n;
	fread(&n, sizeof(int), 1, m_fp);
	if (n>0) fread(sz, sizeof(char), n, m_fp);
	sz[n] = 0;
	return (*this);
}

DumpFile& DumpFile::operator >> (double a[3][3])
{
	fread(a, sizeof(double), 9, m_fp);
	return (*this);
}

DumpFile& DumpFile::operator >> (std::vector<bool>& v)
{
	int n;
	fread(&n, sizeof(int), 1, m_fp);
	if (n > 0) v.resize(n); else v.clear();
	for (int i=0; i<n; ++i)
	{
		int a;
		fread(&a, sizeof(int), 1, m_fp);
		v[i] = (a == 1);
	}
	return (*this);
}
