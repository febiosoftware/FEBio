#include "stdafx.h"
#include "DumpStream.h"

//-----------------------------------------------------------------------------
DumpStream::DumpStream(FEModel& fem) : m_fem(fem)
{
	m_bsave = false;
	m_bshallow = false;
	m_bytes_serialized = 0;
	m_ptr_lock = false;
}

//-----------------------------------------------------------------------------
//! See if the stream is used for input or output
bool DumpStream::IsSaving() const { return m_bsave; }

//-----------------------------------------------------------------------------
//! See if the stream is used for input
bool DumpStream::IsLoading() const { return !m_bsave; }

//-----------------------------------------------------------------------------
//! See if shallow flag is set
bool DumpStream::IsShallow() const { return m_bshallow; }

//-----------------------------------------------------------------------------
DumpStream::~DumpStream()
{
	m_ptr.clear();
	m_bytes_serialized = 0;
}

//-----------------------------------------------------------------------------
void DumpStream::Open(bool bsave, bool bshallow)
{
	m_bsave = bsave;
	m_bshallow = bshallow;
	m_bytes_serialized = 0;
	m_ptr_lock = false;

	// add the "null" pointer
	m_ptr.clear();
	Pointer p = { 0, 0 };
	m_ptr.push_back(p);
}

//-----------------------------------------------------------------------------
void DumpStream::check()
{
	if (IsSaving())
	{
		write(&m_bytes_serialized, sizeof(m_bytes_serialized), 1);
	}
	else
	{
		size_t nsize;
		read(&nsize, sizeof(m_bytes_serialized), 1);
		assert(m_bytes_serialized == nsize);
		if (m_bytes_serialized != nsize) throw DumpStream::ReadError();
	}
}

//-----------------------------------------------------------------------------
void DumpStream::LockPointerTable()
{
	m_ptr_lock = true;
}

//-----------------------------------------------------------------------------
void DumpStream::UnlockPointerTable()
{
	m_ptr_lock = false;
}

//-----------------------------------------------------------------------------
DumpStream& DumpStream::operator << (const char* sz) 
{ 
	int n = (sz ? (int)strlen(sz) : 0);
	m_bytes_serialized += write(&n, sizeof(int), 1);
	if (sz) m_bytes_serialized += write(sz, sizeof(char), n);
	return (*this);
}

//-----------------------------------------------------------------------------
DumpStream& DumpStream::operator << (char* sz) 
{ 
	int n = (sz ? (int)strlen(sz) : 0); 
	m_bytes_serialized += write(&n, sizeof(int), 1);
	if (sz) m_bytes_serialized += write(sz, sizeof(char), n);
	return (*this);
}

//-----------------------------------------------------------------------------
DumpStream& DumpStream::operator<<(std::string& s)
{
	const char* sz = s.c_str();
	this->operator<<(sz);
	return *this;
}

//-----------------------------------------------------------------------------
DumpStream& DumpStream::operator<<(const std::string& s)
{
	const char* sz = s.c_str();
	this->operator<<(sz);
	return *this;
}

//-----------------------------------------------------------------------------
DumpStream& DumpStream::operator << (bool b)
{
	int n = (b ? 1 : 0);
	m_bytes_serialized += write(&n, sizeof(n), 1);
	return *this;
}

//-----------------------------------------------------------------------------
DumpStream& DumpStream::operator << (int n)
{
	m_bytes_serialized += write(&n, sizeof(int), 1);
	return *this;
}

//-----------------------------------------------------------------------------
DumpStream& DumpStream::operator << (const double a[3][3])
{
	m_bytes_serialized += write(a, sizeof(double), 9);
	return (*this);
}

//-----------------------------------------------------------------------------
DumpStream& DumpStream::operator >> (char* sz) 
{ 
	int n;
	m_bytes_serialized += read(&n, sizeof(int), 1);
	if (n>0) m_bytes_serialized += read(sz, sizeof(char), n);
	sz[n] = 0;
	return (*this);
}

//-----------------------------------------------------------------------------
DumpStream& DumpStream::operator >> (std::string& s)
{
	int n;
	m_bytes_serialized += read(&n, sizeof(int), 1);
	char* tmp = new char[n + 1];
	if (n > 0) m_bytes_serialized += read(tmp, sizeof(char), n);
	tmp[n] = 0;
	s = std::string(tmp);
	delete [] tmp;
	return *this;
}

//-----------------------------------------------------------------------------
DumpStream& DumpStream::operator >> (bool& b)
{
	int n;
	m_bytes_serialized += read(&n, sizeof(int), 1);
	b = (n == 1);
	return *this;
}

//-----------------------------------------------------------------------------
DumpStream& DumpStream::operator >> (double a[3][3])
{
	m_bytes_serialized += read(a, sizeof(double), 9);
	return (*this);
}

//-----------------------------------------------------------------------------
int DumpStream::FindPointer(void* p)
{
	for (int i = 0; i < (int)m_ptr.size(); ++i)
	{
		if (m_ptr[i].pd == p) return i;
	}
	return -1;
}

//-----------------------------------------------------------------------------
int DumpStream::FindPointer(int id)
{
	for (int i = 0; i < (int)m_ptr.size(); ++i)
	{
		if (m_ptr[i].id == id) return i;
	}
	return -1;
}

//-----------------------------------------------------------------------------
void DumpStream::AddPointer(void* p)
{
	if (m_ptr_lock) return;
	if (p == nullptr) { assert(false); return;	}
	Pointer ptr;
	ptr.pd = p;
	ptr.id = (int)m_ptr.size();
	m_ptr.push_back(ptr);
}
