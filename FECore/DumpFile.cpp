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
	m_nindex = 0;
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
	m_nindex = 0;
}

//! write buffer to archive
size_t DumpFile::write(const void* pd, size_t size, size_t count)
{
	assert(m_bsave == true);
	m_nindex += (int)(size*count);
	return fwrite(pd, size, count, m_fp);
}

//! read buffer from archive
size_t DumpFile::read(void* pd, size_t size, size_t count)
{
	assert(m_bsave == false);
	m_nindex += (int)(size*count);
	return fread(pd, size, count, m_fp);
}


DumpFile& DumpFile::operator << (const char* sz) 
{ 
	int n = (sz ? strlen(sz) : 0); 
	write(&n, sizeof(int), 1);
	if (sz) write(sz, sizeof(char), n);
	return (*this);
}

DumpFile& DumpFile::operator << (char* sz) 
{ 
	int n = (sz ? strlen(sz) : 0); 
	write(&n, sizeof(int), 1);
	if (sz) write(sz, sizeof(char), n);
	return (*this);
}

DumpFile& DumpFile::operator << (const double a[3][3])
{
	write(a, sizeof(double), 9);
	return (*this);
}

DumpFile& DumpFile::operator << (std::vector<bool>& v)
{
	int n = v.size();
	write(&n, sizeof(int), 1);
	for (int i=0; i<n; ++i) { int a = (v[i]?1:0); write(&a, sizeof(int), 1); }
	return (*this);
}

DumpFile& DumpFile::operator >> (char* sz) 
{ 
	int n;
	read(&n, sizeof(int), 1);
	if (n>0) read(sz, sizeof(char), n);
	sz[n] = 0;
	return (*this);
}

DumpFile& DumpFile::operator >> (double a[3][3])
{
	read(a, sizeof(double), 9);
	return (*this);
}

DumpFile& DumpFile::operator >> (std::vector<bool>& v)
{
	int n;
	read(&n, sizeof(int), 1);
	if (n > 0) v.resize(n); else v.clear();
	for (int i=0; i<n; ++i)
	{
		int a;
		read(&a, sizeof(int), 1);
		v[i] = (a == 1);
	}
	return (*this);
}
