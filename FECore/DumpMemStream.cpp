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
#include "DumpMemStream.h"
#include <assert.h>
#include <memory.h>

//-----------------------------------------------------------------------------
DumpMemStream::DumpMemStream(FEModel& fem) : DumpStream(fem)
{
	m_pb = 0;
	m_pd = 0;
	m_nsize = 0;
	m_nreserved = 0;

	Open(true, true);
}

//-----------------------------------------------------------------------------
void DumpMemStream::clear()
{
	delete [] m_pb;
	m_pb = 0;
	m_pd = 0;
	m_nsize = 0;
	m_nreserved = 0;

	// Since we can't read from an empty stream
	// we restore write mode.
	Open(true, true);
}

//-----------------------------------------------------------------------------
void DumpMemStream::Open(bool bsave, bool bshallow)
{
	DumpStream::Open(bsave, bshallow);
	if (m_pb) set_position(0);
}

//-----------------------------------------------------------------------------
bool DumpMemStream::EndOfStream() const
{
	return (bytesSerialized() >= m_nsize);
}

//-----------------------------------------------------------------------------
DumpMemStream::~DumpMemStream()
{
	clear();
}

//-----------------------------------------------------------------------------
void DumpMemStream::set_position(size_t l)
{
	assert((l >= 0) && (l < m_nreserved));
	m_pd = m_pb + l;
}

//-----------------------------------------------------------------------------
void DumpMemStream::grow_buffer(size_t l)
{
	if (l <= 0) return;

	char* pnew = new char[m_nreserved + l];
	if (m_pb)
	{
		memcpy(pnew, m_pb, m_nreserved);
		delete [] m_pb;
	}
	m_pb = pnew;
	m_pd = m_pb + m_nsize;
	m_nreserved += l;
}

//-----------------------------------------------------------------------------
size_t DumpMemStream::write(const void* pd, size_t size, size_t count)
{
	assert(IsSaving());
	size_t nsize = count*size;
	size_t lpos = (size_t)(m_pd - m_pb);
	if (lpos + nsize > m_nreserved) grow_buffer(nsize + 3*m_nreserved/2);
	memcpy(m_pd, pd, nsize);

	m_pd += nsize;
	lpos += nsize;
	if (lpos > m_nsize) m_nsize = lpos;

	return nsize;
}

//-----------------------------------------------------------------------------
size_t DumpMemStream::read(void* pd, size_t size, size_t count)
{
	assert(IsSaving()==false);
	size_t nsize = count*size;
	memcpy(pd, m_pd, nsize);
	m_pd += nsize;
	return nsize;
}
