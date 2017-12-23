#include "stdafx.h"
#include "DumpStream.h"

//-----------------------------------------------------------------------------
DumpStream::DumpStream(FEModel& fem) : m_fem(fem)
{
	m_bsave = false;
	m_bshallow = false;
}

//-----------------------------------------------------------------------------
DumpStream::~DumpStream()
{
}

//-----------------------------------------------------------------------------
void DumpStream::Open(bool bsave, bool bshallow)
{
	m_bsave = bsave;
	m_bshallow = bshallow;
}

//-----------------------------------------------------------------------------
DumpStream& DumpStream::operator << (const char* sz) 
{ 
	int n = (sz ? (int)strlen(sz) : 0);
	write(&n, sizeof(int), 1);
	if (sz) write(sz, sizeof(char), n);
	return (*this);
}

//-----------------------------------------------------------------------------
DumpStream& DumpStream::operator << (char* sz) 
{ 
	int n = (sz ? (int)strlen(sz) : 0); 
	write(&n, sizeof(int), 1);
	if (sz) write(sz, sizeof(char), n);
	return (*this);
}

//-----------------------------------------------------------------------------
DumpStream& DumpStream::operator<<(const std::string& s)
{
	const char* sz = s.c_str();
	this->operator<<(sz);
	return *this;
}

//-----------------------------------------------------------------------------
DumpStream& DumpStream::operator << (const double a[3][3])
{
	write(a, sizeof(double), 9);
	return (*this);
}

//-----------------------------------------------------------------------------
DumpStream& DumpStream::operator >> (char* sz) 
{ 
	int n;
	read(&n, sizeof(int), 1);
	if (n>0) read(sz, sizeof(char), n);
	sz[n] = 0;
	return (*this);
}

//-----------------------------------------------------------------------------
DumpStream& DumpStream::operator >> (std::string& s)
{
	char buf[64]={0};
	this->operator>>(buf);
	s = std::string(buf);
	return *this;
}

//-----------------------------------------------------------------------------
DumpStream& DumpStream::operator >> (double a[3][3])
{
	read(a, sizeof(double), 9);
	return (*this);
}
