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
DumpMemStream::~DumpMemStream()
{
	clear();
}

//-----------------------------------------------------------------------------
void DumpMemStream::set_position(int l)
{
	assert((l >= 0) && (l < m_nreserved));
	m_pd = m_pb + l;
	m_nsize = l;
}

//-----------------------------------------------------------------------------
void DumpMemStream::grow_buffer(int l)
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
	int nsize = (int)(count*size);
	if (m_nsize + nsize > m_nreserved) grow_buffer(nsize + m_nreserved/10);
	memcpy(m_pd, pd, nsize);
	m_pd += nsize;
	m_nsize += nsize;
	return count;
}

//-----------------------------------------------------------------------------
size_t DumpMemStream::read(void* pd, size_t size, size_t count)
{
	assert(IsSaving()==false);
	int nsize = (int)(count*size);
	memcpy(pd, m_pd, nsize);
	m_pd += nsize;
	m_nsize += nsize;
	return count;
}

//-----------------------------------------------------------------------------
void DumpMemStream::check()
{
	if (IsSaving())
	{
		write(&m_nsize, sizeof(m_nsize), 1);
	}
	else
	{
		int n = m_nsize, nsize;
		read(&nsize, sizeof(m_nsize), 1);
		assert(n==nsize);
		if (n != nsize) throw DumpStream::ReadError();
	}
}
