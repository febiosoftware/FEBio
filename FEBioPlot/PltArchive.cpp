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
#include "PltArchive.h"
#include <assert.h>

#ifdef HAVE_ZLIB
#include "zlib.h"
static z_stream strm;
#endif

//=============================================================================
// FileStream
//=============================================================================
FileStream::FileStream(FILE* fp, bool owner)
{
	m_bufsize = 262144;	// = 256K
	m_current = 0;
	m_buf  = new unsigned char[m_bufsize];
	m_pout = new unsigned char[m_bufsize];
	m_ncompress = 0;
	m_fp = fp;
	m_fileOwner = owner;
}

FileStream::~FileStream()
{
	Close();
	delete [] m_buf;
	delete [] m_pout;
	m_buf = 0;
	m_pout = 0;
}

bool FileStream::Open(const char* szfile)
{
	m_fp = fopen(szfile, "rb");
	if (m_fp == 0) return false;
	return true;
}

bool FileStream::Append(const char* szfile)
{
	m_fp = fopen(szfile, "a+b");
	return (m_fp != 0);
}

bool FileStream::Create(const char* szfile)
{
	m_fp = fopen(szfile, "wb");
	return (m_fp != 0);
}

void FileStream::Close()
{
	if (m_fp)
	{
		Flush();
		if (m_fileOwner) fclose(m_fp);
	}
	m_fp = 0;
}

void FileStream::BeginStreaming()
{
#ifdef HAVE_ZLIB
	if (m_ncompress)
	{
		strm.zalloc = Z_NULL;
		strm.zfree = Z_NULL;
		strm.opaque = Z_NULL;
		deflateInit(&strm, -1);
	}
#endif
}

void FileStream::EndStreaming()
{
	Flush();
#ifdef HAVE_ZLIB
	if (m_ncompress)
	{
		strm.avail_in = 0;
		strm.next_in = 0;

		/* run deflate() on input until output buffer not full, finish
		compression if all of source has been read in */
		do {
			strm.avail_out = m_bufsize;
			strm.next_out = m_pout;
			int ret = deflate(&strm, Z_FINISH);    /* no bad return value */
			assert(ret != Z_STREAM_ERROR);  /* state not clobbered */
			int have = m_bufsize - strm.avail_out;
			fwrite(m_pout, 1, have, m_fp);
		} while (strm.avail_out == 0);
		assert(strm.avail_in == 0);     /* all input will be used */

		// all done
		deflateEnd(&strm);

		fflush(m_fp);
	}
#endif
}

void FileStream::Write(void* pd, size_t Size, size_t Count)
{
	unsigned char* pdata = (unsigned char*) pd;
	size_t nsize = Size*Count;
	while (nsize > 0)
	{
		if (m_current + nsize < m_bufsize)
		{
			memcpy(m_buf + m_current, pdata, nsize);
			m_current += nsize;
			nsize = 0;
		}
		else
		{
			int nblock = m_bufsize - m_current;
			if (nblock>0) { memcpy(m_buf + m_current, pdata, nblock); m_current += nblock; }
			Flush();
			pdata += nblock;
			nsize -= nblock;
		}
	}
}

void FileStream::Flush()
{
#ifdef HAVE_ZLIB
	if (m_ncompress)
	{
		strm.avail_in = m_current;
		strm.next_in = m_buf;

		/* run deflate() on input until output buffer not full, finish
		compression if all of source has been read in */
		do {
			strm.avail_out = m_bufsize;
			strm.next_out = m_pout;
			int ret = deflate(&strm, Z_NO_FLUSH);    /* no bad return value */
			assert(ret != Z_STREAM_ERROR);  /* state not clobbered */
			int have = m_bufsize - strm.avail_out;
			fwrite(m_pout, 1, have, m_fp);
		} while (strm.avail_out == 0);
		assert(strm.avail_in == 0);     /* all input will be used */
	}
	else
	{
		if (m_fp) fwrite(m_buf, m_current, 1, m_fp);
	}
#else
	if (m_fp) fwrite(m_buf, m_current, 1, m_fp);
#endif

	// flush the file
	if (m_fp) fflush(m_fp);

	// reset current data pointer
	m_current = 0;
}

size_t FileStream::read(void* pd, size_t Size, size_t Count)
{
	return fread(pd, Size, Count, m_fp);
}

long FileStream::tell()
{
	return ftell(m_fp);
}

void FileStream::seek(long noff, int norigin)
{
	fseek(m_fp, noff, norigin);
}


//=============================================================================
// PltArchive
//=============================================================================

PltArchive::PltArchive()
{
	m_fp = 0;
	m_pRoot = 0;
	m_pChunk = 0;
	m_bSaving = true;
}

PltArchive::~PltArchive()
{
	Close();
}

void PltArchive::Close()
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
	if (m_fp)
	{
		m_fp->Close();
		delete m_fp;
		m_fp = 0;
	}
}

void PltArchive::SetCompression(int n)
{
	if (m_fp) m_fp->SetCompression(n);
}

void PltArchive::Flush()
{
	if (m_fp && m_pRoot)
	{
		m_fp->BeginStreaming();
		m_pRoot->Write(m_fp);
		m_fp->EndStreaming();
	}
	delete m_pRoot;
	m_pRoot = 0;
	m_pChunk = 0;
}

bool PltArchive::Create(const char* szfile)
{
	// attempt to create the file
	assert(m_fp == 0);
	m_fp = new FileStream();
	if (m_fp->Create(szfile) == false) return false;

	// write the root tag 
	unsigned int ntag = 0x00464542;
	m_fp->Write(&ntag, sizeof(int), 1);

	m_bSaving = true;

	return true;
}

void PltArchive::BeginChunk(unsigned int id)
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

void PltArchive::EndChunk()
{
	if (m_pChunk != m_pRoot)
		m_pChunk = m_pChunk->GetParent();
	else 
	{
		Flush();
	}
}



//-----------------------------------------------------------------------------

bool PltArchive::Open(const char* szfile)
{
	// try to open the file
	assert(m_fp == 0);
	m_fp = new FileStream();
	if (m_fp->Open(szfile) == false) return false;

	// read the root tag
	unsigned int ntag;
	m_fp->read(&ntag, sizeof(int), 1);
	if (ntag != 0x00464542) { Close(); return false; }

	m_bSaving = false;
	m_bend = false;
	
	return true;
}

bool PltArchive::Append(const char* szfile)
{
	// reopen the plot file for appending
	assert(m_fp == 0);
	m_fp = new FileStream();
	if (m_fp->Append(szfile) == false) return false;
	m_bSaving = true;
	return true;
}

int PltArchive::OpenChunk()
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
	pc->lpos = m_fp->tell();

	// add it to the stack
	m_Chunk.push(pc);

	return IO_OK;
}

void PltArchive::CloseChunk()
{
	// pop the last chunk
	CHUNK* pc = m_Chunk.top(); m_Chunk.pop();

	// get the current file position
	long lpos = m_fp->tell();

	// calculate the offset to the end of the chunk
	int noff = pc->nsize - (lpos - pc->lpos);

	// skip any remaining part in the chunk
	// I wonder if this can really happen
	if (noff != 0)
	{
		m_fp->seek(noff, SEEK_CUR);
		lpos = m_fp->tell();
	}

	// delete this chunk
	delete pc;

	// take a peek at the parent
	if (m_Chunk.empty())
	{
		// we just deleted the root chunk
		m_bend = true;
	}
	else
	{
		pc = m_Chunk.top();
		int noff = pc->nsize - (lpos - pc->lpos);
		if (noff == 0) m_bend = true;
	}
}

unsigned int PltArchive::GetChunkID()
{
	CHUNK* pc = m_Chunk.top();
	assert(pc);
	return pc->id;
}
