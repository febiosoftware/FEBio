// Archive.cpp: implementation of the Archive class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "Archive.h"

//////////////////////////////////////////////////////////////////////
// Archive
//////////////////////////////////////////////////////////////////////

Archive::Archive()
{
	m_fp = 0;
	m_pRoot = 0;
	m_pChunk = 0;
}

Archive::~Archive()
{
	Close();
}

void Archive::Close()
{
	if (m_pRoot) Flush();
	if (m_fp) fclose(m_fp);
	m_fp = 0;
}

void Archive::Flush()
{
	if (m_fp && m_pRoot) m_pRoot->Write(m_fp);
	delete m_pRoot;
	m_pRoot = 0;
	m_pChunk = 0;
}

bool Archive::Create(const char* szfile)
{
	// attempt to create the file
	m_fp = fopen(szfile, "wb");
	if (m_fp == 0) return false;

	// write the master tag 
	unsigned int ntag = 0x00464542;
	fwrite(&ntag, sizeof(int), 1, m_fp);

	return true;
}

void Archive::BeginChunk(unsigned int id)
{
	if (m_pRoot == 0)
	{
		m_pRoot = new OBranch(id);
		m_pChunk = m_pRoot;
	}
	else
	{
		// create a new branch
		OBranch* pbranch = new OBranch(id);

		// attach it to the current branch
		m_pChunk->AddChild(pbranch);

		// move the current branch pointer
		m_pChunk = pbranch;
	}
}

void Archive::EndChunk()
{
	if (m_pChunk != m_pRoot)
		m_pChunk = m_pChunk->GetParent();
	else 
	{
		Flush();
	}
}
