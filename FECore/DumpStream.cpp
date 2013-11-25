#include "stdafx.h"
#include "DumpStream.h"
#include <assert.h>
#include <memory.h>

//-----------------------------------------------------------------------------
DumpStream::DumpStream()
{
	m_pb = 0;
	m_pd = 0;
	m_nsize = 0;
	m_nreserved = 0;
}

//-----------------------------------------------------------------------------
void DumpStream::clear()
{
	delete [] m_pb;
	m_pb = 0;
	m_pd = 0;
	m_nsize = 0;
	m_nreserved = 0;
}

//-----------------------------------------------------------------------------
DumpStream::~DumpStream()
{
	clear();
}

//-----------------------------------------------------------------------------
void DumpStream::set_position(int l)
{
	assert((l >= 0) && (l < m_nreserved));
	m_pd = m_pb + l;
	m_nsize = l;
}

//-----------------------------------------------------------------------------
void DumpStream::grow_buffer(int l)
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
void DumpStream::write(void* pd, int nsize)
{
	if (m_nsize + nsize > m_nreserved) grow_buffer(nsize + m_nreserved/10);
	memcpy(m_pd, pd, nsize);
	m_pd += nsize;
	m_nsize += nsize;
}

//-----------------------------------------------------------------------------
void DumpStream::read(void* pd, int nsize)
{
	memcpy(pd, m_pd, nsize);
	m_pd += nsize;
	m_nsize += nsize;
}
