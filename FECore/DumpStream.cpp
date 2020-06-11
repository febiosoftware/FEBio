/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2020 University of Utah, The Trustees of Columbia University in 
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
#include "DumpStream.h"
#include "matrix.h"

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
	assert(FindPointer(p) == -1);
	Pointer ptr;
	ptr.pd = p;
	ptr.id = (int)m_ptr.size();
	m_ptr.push_back(ptr);
}

//-----------------------------------------------------------------------------
DumpStream& DumpStream::write_matrix(matrix& o)
{
	DumpStream& ar = *this;
	int nr = o.rows();
	int nc = o.columns();
	ar << nr << nc;
	int nsize = nr*nc;
	if (nsize > 0)
	{
		vector<double> data;
		data.reserve(nr*nc);
		for (int i = 0; i < nr; ++i)
			for (int j = 0; j < nc; ++j) data.push_back(o(i, j));

		ar << data;
	}
	return *this;
}

//-----------------------------------------------------------------------------
DumpStream& DumpStream::read_matrix(matrix& o)
{
	DumpStream& ar = *this;
	int nr = 0, nc = 0;
	ar >> nr >> nc;
	int nsize = nr*nc;
	if (nsize > 0)
	{
		o.resize(nr, nc);
		vector<double> data;
		ar >> data;
		int n = 0;
		for (int i = 0; i < nr; ++i)
			for (int j = 0; j < nc; ++j) o(i, j) = data[n++];
	}

	return *this;
}
