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
#include "DumpStream.h"
#include "matrix.h"

//-----------------------------------------------------------------------------
DumpStream::DumpStream(FEModel& fem) : m_fem(fem)
{
	m_bsave = false;
	m_bshallow = false;
	m_bytes_serialized = 0;
	m_ptr_lock = false;

#ifdef _DEBUG
	m_btypeInfo = false;
#else
	m_btypeInfo = false;
#endif
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
	m_ptrOut.clear();
	m_ptrIn.clear();
	m_bytes_serialized = 0;
}

//-----------------------------------------------------------------------------
// set the write type info flag
void DumpStream::WriteTypeInfo(bool b)
{
	m_btypeInfo = b;
}

//-----------------------------------------------------------------------------
// see if the stream has type info
bool DumpStream::HasTypeInfo() const
{
	return m_btypeInfo;
}

//-----------------------------------------------------------------------------
void DumpStream::Open(bool bsave, bool bshallow)
{
	m_bsave = bsave;
	m_bshallow = bshallow;
	m_bytes_serialized = 0;
	m_ptr_lock = false;

	// add the "null" pointer
	if (bsave)
	{
		m_ptrOut.clear();
		m_ptrOut[nullptr] = 0;
	}
	else
	{
		m_ptrIn.clear();
		m_ptrIn.push_back(nullptr);
	}
}

//-----------------------------------------------------------------------------
void DumpStream::check()
{
	if (IsSaving())
	{
		m_bytes_serialized += write(&m_bytes_serialized, sizeof(m_bytes_serialized), 1);
	}
	else
	{
		size_t nsize;
		size_t tmp = read(&nsize, sizeof(m_bytes_serialized), 1);
		assert(m_bytes_serialized == nsize);
		if (m_bytes_serialized != nsize) throw DumpStream::ReadError();
		m_bytes_serialized += tmp;
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
	if (sz && (n > 0)) m_bytes_serialized += write(sz, sizeof(char), n);
	return (*this);
}

//-----------------------------------------------------------------------------
DumpStream& DumpStream::operator << (char* sz) 
{ 
	int n = (sz ? (int)strlen(sz) : 0); 
	m_bytes_serialized += write(&n, sizeof(int), 1);
	if (sz && (n > 0)) m_bytes_serialized += write(sz, sizeof(char), n);
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
	if (m_btypeInfo) writeType(TypeID::TYPE_INT);
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
	assert(IsSaving());
	auto it = m_ptrOut.find(p);
	if (it != m_ptrOut.end()) return it->second;
	return -1;
}

//-----------------------------------------------------------------------------
void DumpStream::AddPointer(void* p)
{
	if (m_ptr_lock) return;
	if (p == nullptr) { assert(false); return;	}

	if (IsSaving())
	{
		assert(FindPointer(p) == -1);
		int id = (int)m_ptrOut.size();
		m_ptrOut[p] = id;
	}
	else
	{
		m_ptrIn.push_back(p);
	}
}

//-----------------------------------------------------------------------------
DumpStream& DumpStream::write_matrix(matrix& o)
{
	if (m_btypeInfo) writeType(TypeID::TYPE_MATRIX);

	// don't write type info for all components
	bool oldTypeFlag = m_btypeInfo;
	m_btypeInfo = false;

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

	m_btypeInfo = oldTypeFlag;

	return *this;
}

//-----------------------------------------------------------------------------
DumpStream& DumpStream::read_matrix(matrix& o)
{
	if (m_btypeInfo) readType(TypeID::TYPE_MATRIX);

	// don't read type info for all components
	bool oldTypeFlag = m_btypeInfo;
	m_btypeInfo = false;

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

	m_btypeInfo = oldTypeFlag;

	return *this;
}

//-----------------------------------------------------------------------------
// read the next block
bool DumpStream::readBlock(DataBlock& d)
{
	// make sure we have type info
	if (m_btypeInfo == false) return false;

	// see if we have reached the end of the stream
	if (EndOfStream()) return false;

	// read the data type
	d.m_type = readType();

	// turn off type flag since we already read it
	m_btypeInfo = false;

	// read/allocate data
	switch (d.m_type)
	{
	case TypeID::TYPE_INT     : { int          v; read_raw(v); d.m_pd = new int         (v); } break;
	case TypeID::TYPE_UINT    : { unsigned int v; read_raw(v); d.m_pd = new unsigned int(v); } break;
	case TypeID::TYPE_FLOAT   : { float        v; read_raw(v); d.m_pd = new float       (v); } break;
	case TypeID::TYPE_DOUBLE  : { double       v; read_raw(v); d.m_pd = new double      (v); } break;
	case TypeID::TYPE_VEC2D   : { vec2d        v; read_raw(v); d.m_pd = new vec2d       (v); } break;
	case TypeID::TYPE_VEC3D   : { vec3d        v; read_raw(v); d.m_pd = new vec3d       (v); } break;
	case TypeID::TYPE_MAT2D   : { mat2d        v; read_raw(v); d.m_pd = new mat2d       (v); } break;
	case TypeID::TYPE_MAT3D   : { mat3d        v; read_raw(v); d.m_pd = new mat3d       (v); } break;
	case TypeID::TYPE_MAT3DD  : { mat3dd       v; read_raw(v); d.m_pd = new mat3dd      (v); } break;
	case TypeID::TYPE_MAT3DS  : { mat3ds       v; read_raw(v); d.m_pd = new mat3ds      (v); } break;
	case TypeID::TYPE_MAT3DA  : { mat3da       v; read_raw(v); d.m_pd = new mat3da      (v); } break;
	case TypeID::TYPE_QUATD   : { quatd        v; read_raw(v); d.m_pd = new quatd       (v); } break;
	case TypeID::TYPE_TENS3DS : { tens3ds      v; read_raw(v); d.m_pd = new tens3ds     (v); } break;
	case TypeID::TYPE_TENS3DRS: { tens3drs     v; read_raw(v); d.m_pd = new tens3drs    (v); } break;
	case TypeID::TYPE_MATRIX  : { matrix       v; read_raw(v); d.m_pd = new matrix      (v); } break;
	default:
		assert(false);
		m_btypeInfo = true;
		return false;
	}

	// turn the type info flag back on
	m_btypeInfo = true;
	return true;
}
