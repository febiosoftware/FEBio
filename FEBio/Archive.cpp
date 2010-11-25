// Archive.cpp: implementation of the Archive class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "Archive.h"


Archive::Archive()
{
	m_fp = 0;
}

Archive::~Archive()
{
	Close();
}

void Archive::Close()
{
	if (m_fp) fclose(m_fp);
	m_fp = 0;
}

bool Archive::Create(const char* sz)
{
	if (m_fp) Close();

	m_fp = fopen(sz, "wb");
	if (m_fp == 0) return false;

	return true;
}

void Archive::BeginChunk(unsigned int id, int nsize)
{
	// write the chunk ID
	unsigned int n[2] = {id, nsize};
	fwrite(n, sizeof(unsigned int), 2, m_fp);
	m_tag.push(id);
}

void Archive::EndChunk()
{
	assert(m_tag.empty() == false);

	unsigned int nid = m_tag.front();
	m_tag.pop();
	unsigned int n[2] = {nid, (unsigned int) (-1)};
	fwrite(n, sizeof(unsigned int), 2, m_fp);
}
