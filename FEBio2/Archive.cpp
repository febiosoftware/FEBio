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
	m_bSaving = true;
}

Archive::~Archive()
{
	Close();
}

void Archive::Close()
{
	if (m_bSaving)
	{
		if (m_pRoot) Flush();
	}
	else 
	{
		while (m_Chunk.empty() == false) CloseChunk();
		m_bend = true;
	}

	// close the file
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

	m_bSaving = true;

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



//-----------------------------------------------------------------------------

bool Archive::Open(const char* szfile)
{
	// try to open the file
	m_fp = fopen(szfile, "rb");
	if (m_fp == 0) return false;

	// read the master tag
	unsigned int ntag;
	fread(&ntag, sizeof(int), 1, m_fp);
	if (ntag != 0x00464542) { Close(); return false; }

	m_bSaving = false;
	m_bend = false;
	
	return true;
}

bool Archive::Append(const char* szfile)
{
	// reopen the plot file for appending
	m_fp = fopen(szfile, "a+b");
	return (m_fp != 0);
}

int Archive::OpenChunk()
{
	// see if the end flag was set
	// in that case we first need to clear the flag
	if (m_bend)
	{
		m_bend = false;
		return IO_END;
	}

	// create a new chunk
	CHUNK* pc = new CHUNK;

	// read the chunk ID
	read(pc->id);

	// read the chunk size
	read(pc->nsize);

	if (pc->nsize == 0) m_bend = true;

	// record the position
	pc->lpos = ftell(m_fp);

	// add it to the stack
	m_Chunk.push(pc);

	return IO_OK;
}

void Archive::CloseChunk()
{
	// pop the last chunk
	CHUNK* pc = m_Chunk.top(); m_Chunk.pop();

	// get the current file position
	long lpos = ftell(m_fp);

	// calculate the offset to the end of the chunk
	int noff = pc->nsize - (lpos - pc->lpos);

	// skip any remaining part in the chunk
	// I wonder if this can really happen
	if (noff != 0)
	{
		fseek(m_fp, noff, SEEK_CUR);
		lpos = ftell(m_fp);
	}

	// delete this chunk
	delete pc;

	// take a peek at the parent
	if (m_Chunk.empty())
	{
		// we just deleted the master chunk
		m_bend = true;
	}
	else
	{
		pc = m_Chunk.top();
		int noff = pc->nsize - (lpos - pc->lpos);
		if (noff == 0) m_bend = true;
	}
}

unsigned int Archive::GetChunkID()
{
	CHUNK* pc = m_Chunk.top();
	assert(pc);
	return pc->id;
}
