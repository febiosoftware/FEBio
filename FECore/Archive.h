// Archive.h: interface for the Archive class.
//
//////////////////////////////////////////////////////////////////////
#pragma once
#include <assert.h>
#include <string.h>
#include <list>
#include <vector>
#include <stack>
using namespace std;

//-----------------------------------------------------------------------------
//! helper class for writing buffered data to file
class FileStream
{
public:
	FileStream();
	~FileStream();

	bool Create(const char* szfile);
	bool Open(const char* szfile);
	bool Append(const char* szfile);
	void Close();

	void Write(void* pd, size_t Size, size_t Count);

	void Flush();

	// \todo temporary reading functions. Needs to be replaced with buffered functions
	size_t read(void* pd, size_t Size, size_t Count);
	long tell();
	void seek(long noff, int norigin);

	void BeginStreaming();
	void EndStreaming();

	void SetCompression(int n) { m_ncompress = n; }

private:
	FILE*	m_fp;
	size_t	m_bufsize;		//!< buffer size
	size_t	m_current;		//!< current index
	unsigned char*	m_buf;	//!< buffer
	unsigned char*	m_pout;	//!< temp buffer when writing
	int		m_ncompress;	//!< compression level
};


//----------------------
// Output archive

class OBranch;

class OChunk
{
public:
	OChunk(unsigned int nid) { m_nID = nid; m_pParent = 0; }
	virtual ~OChunk(){}

	unsigned int GetID() { return m_nID; }

	virtual void Write(FileStream* fp) = 0;
	virtual int Size() = 0;

	void SetParent(OBranch* pparent) { m_pParent = pparent; }
	OBranch* GetParent() { return m_pParent; }

protected:
	int			m_nID;
	OBranch*	m_pParent;
};

class OBranch : public OChunk
{
public:
	OBranch(unsigned int nid) : OChunk(nid) {}
	~OBranch()
	{
		list<OChunk*>::iterator pc;
		for (pc = m_child.begin(); pc != m_child.end(); ++pc) delete (*pc);
		m_child.clear();
	}

	int Size()
	{
		int nsize = 0;
		list<OChunk*>::iterator pc;
		for (pc = m_child.begin(); pc != m_child.end(); ++pc) nsize += (*pc)->Size() + 2*sizeof(unsigned int);
		return nsize;
	}

	void Write(FileStream* fp)
	{
		fp->Write(&m_nID  , sizeof(unsigned int), 1);

		unsigned int nsize = Size();
		fp->Write(&nsize, sizeof(unsigned int), 1);

		list<OChunk*>::iterator pc;
		for (pc = m_child.begin(); pc != m_child.end(); ++pc) (*pc)->Write(fp);
	}

	void AddChild(OChunk* pc) { m_child.push_back(pc); pc->SetParent(this); }

protected:
	list<OChunk*>	m_child;
};

template <typename T>
class OLeaf : public OChunk
{
public:
	OLeaf(unsigned int nid, const T& d) : OChunk(nid) { m_d = d; }

	int Size() { return sizeof(T); }

	void Write(FileStream* fp)
	{
		fp->Write(&m_nID  , sizeof(unsigned int), 1);
		unsigned int nsize = sizeof(T);
		fp->Write(&nsize, sizeof(unsigned int), 1);
		fp->Write(&m_d, sizeof(T), 1);
	}

protected:
	T	m_d;
};

template <typename T>
class OLeaf<T*> : public OChunk
{
public:
	OLeaf(unsigned int nid, const T* pd, int nsize) : OChunk(nid)
	{
		assert(nsize > 0);
		m_pd = new T[nsize];
		memcpy(m_pd, pd, sizeof(T)*nsize);
		m_nsize = nsize;
	}
	~OLeaf() { delete m_pd; }

	int Size() { return sizeof(T)*m_nsize; }
	void Write(FileStream* fp)
	{
		fp->Write(&m_nID , sizeof(unsigned int), 1);
		unsigned int nsize = Size();
		fp->Write(&nsize , sizeof(unsigned int), 1);
		fp->Write(m_pd   , sizeof(T), m_nsize);
	}

protected:
	T*		m_pd;
	int		m_nsize;
};

template <>
class OLeaf<const char*> : public OChunk
{
public:
	OLeaf(unsigned int nid, const char* sz) : OChunk(nid)
	{
		int l = strlen(sz);
		m_psz = new char[l+1];
		memcpy(m_psz, sz, l+1);
	}
	~OLeaf() { delete m_psz; }

	int Size() { return strlen(m_psz)+sizeof(int); }
	void Write(FileStream* fp)
	{
		fp->Write(&m_nID , sizeof(unsigned int), 1);
		unsigned int nsize = Size();
		fp->Write(&nsize , sizeof(unsigned int), 1);
		int l = nsize - sizeof(int);
		fp->Write(&l, sizeof(int), 1);
		fp->Write(m_psz, sizeof(char), l);
	}

protected:
	char*	m_psz;
};

template <typename T>
class OLeaf<vector<T> > : public OChunk
{
public:
	OLeaf(unsigned int nid, vector<T>& a) : OChunk(nid)
	{
		m_nsize = a.size();
		assert(m_nsize > 0);
		m_pd = new T[m_nsize];
		memcpy(m_pd, &a[0], sizeof(T)*m_nsize);
	}
	~OLeaf() { delete m_pd; }

	int Size() { return sizeof(T)*m_nsize; }
	void Write(FileStream* fp)
	{
		fp->Write(&m_nID , sizeof(unsigned int), 1);
		unsigned int nsize = Size();
		fp->Write(&nsize , sizeof(unsigned int), 1);
		fp->Write(m_pd   , sizeof(T), m_nsize);
	}

protected:
	T*		m_pd;
	int		m_nsize;
};

enum IOResult { IO_ERROR, IO_OK, IO_END };

class Archive  
{
protected:
	// CHUNK data structure for reading
	struct CHUNK
	{
		unsigned int	id;		// chunk ID
		unsigned int	lpos;	// file position
		unsigned int	nsize;	// size of chunk
	};

public:
	Archive();
	virtual ~Archive();

	// Close archive
	void Close();

	// flush data to file
	void Flush();

	// --- Writing ---

	// Open for writing
	bool Create(const char* szfile);

	// begin a chunk
	void BeginChunk(unsigned int id);

	// end a chunck
	void EndChunk();

	template <typename T> void WriteChunk(unsigned int nid, T& o)
	{
		m_pChunk->AddChild(new OLeaf<T>(nid, o));
	}

	void WriteChunk(unsigned int nid, char* sz)
	{
		m_pChunk->AddChild(new OLeaf<const char*>(nid, sz));
	}

	template <typename T> void WriteChunk(unsigned int nid, T* po, int n)
	{
		m_pChunk->AddChild(new OLeaf<T*>(nid, po, n));
	}

	template <typename T> void WriteChunk(unsigned int nid, vector<T>& a)
	{
		m_pChunk->AddChild(new OLeaf<vector<T> >(nid, a));
	}

	// --- Reading ---

	// Open for reading
	bool Open(const char* sfile);
	bool Append(const char* szfile);

	// Open a chunk
	int OpenChunk();

	// Get the current chunk ID
	unsigned int GetChunkID();

	// Close a chunk
	void CloseChunk();

	// input functions
	IOResult read(char&   c) { size_t nr = m_fp->read(&c, sizeof(char  ), 1); if (nr != 1) return IO_ERROR; return IO_OK; }
	IOResult read(int&    n) { size_t nr = m_fp->read(&n, sizeof(int   ), 1); if (nr != 1) return IO_ERROR; return IO_OK; }
	IOResult read(bool&   b) { size_t nr = m_fp->read(&b, sizeof(bool  ), 1); if (nr != 1) return IO_ERROR; return IO_OK; }
	IOResult read(float&  f) { size_t nr = m_fp->read(&f, sizeof(float ), 1); if (nr != 1) return IO_ERROR; return IO_OK; }
	IOResult read(double& g) { size_t nr = m_fp->read(&g, sizeof(double), 1); if (nr != 1) return IO_ERROR; return IO_OK; }

	IOResult read(unsigned int& n) { size_t nr = m_fp->read(&n, sizeof(unsigned int), 1); if (nr != 1) return IO_ERROR; return IO_OK; }

	IOResult read(char*   pc, int n) { size_t nr = m_fp->read(pc, sizeof(char  ), n); if (nr != n) return IO_ERROR; return IO_OK; }
	IOResult read(int*    pi, int n) { size_t nr = m_fp->read(pi, sizeof(int   ), n); if (nr != n) return IO_ERROR; return IO_OK; }
	IOResult read(bool*   pb, int n) { size_t nr = m_fp->read(pb, sizeof(bool  ), n); if (nr != n) return IO_ERROR; return IO_OK; }
	IOResult read(float*  pf, int n) { size_t nr = m_fp->read(pf, sizeof(float ), n); if (nr != n) return IO_ERROR; return IO_OK; }
	IOResult read(double* pg, int n) { size_t nr = m_fp->read(pg, sizeof(double), n); if (nr != n) return IO_ERROR; return IO_OK; }

	IOResult read(char* sz)
	{
		IOResult ret;
		int l;
		ret = read(l); if (ret != IO_OK) return ret;
		size_t nr = m_fp->read(sz, 1, l); if (nr != l) return IO_ERROR;
		sz[l] = 0;
		return IO_OK;
	}

	void SetCompression(int n);

protected:
	FileStream*	m_fp;		// pointer to file stream
	bool		m_bSaving;	// read or write mode?

	// write data
	OBranch*	m_pRoot;	// chunk tree root
	OBranch*	m_pChunk;	// current chunk

	// read data
	bool			m_bend;		// chunk end flag
	stack<CHUNK*>	m_Chunk;
};
