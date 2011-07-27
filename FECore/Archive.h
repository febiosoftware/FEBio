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

//----------------------
// Output archive

class OBranch;

class OChunk
{
public:
	OChunk(unsigned int nid) { m_nID = nid; m_pParent = 0; }
	virtual ~OChunk(){}

	unsigned int GetID() { return m_nID; }

	virtual void Write(FILE* fp) = 0;
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

	void Write(FILE* fp)
	{
		fwrite(&m_nID  , sizeof(unsigned int), 1, fp);

		unsigned int nsize = Size();
		fwrite(&nsize, sizeof(unsigned int), 1, fp);

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

	void Write(FILE* fp)
	{
		fwrite(&m_nID  , sizeof(unsigned int), 1, fp);
		unsigned int nsize = sizeof(T);
		fwrite(&nsize, sizeof(unsigned int), 1, fp);
		fwrite(&m_d, sizeof(T), 1, fp);
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
	void Write(FILE* fp)
	{
		fwrite(&m_nID , sizeof(unsigned int), 1, fp);
		unsigned int nsize = Size();
		fwrite(&nsize , sizeof(unsigned int), 1, fp);
		fwrite(m_pd   , sizeof(T), m_nsize, fp);
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
	void Write(FILE* fp)
	{
		fwrite(&m_nID , sizeof(unsigned int), 1, fp);
		unsigned int nsize = Size();
		fwrite(&nsize , sizeof(unsigned int), 1, fp);
		int l = nsize - sizeof(int);
		fwrite(&l, sizeof(int), 1, fp);
		fwrite(m_psz, sizeof(char), l, fp);
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
	void Write(FILE* fp)
	{
		fwrite(&m_nID , sizeof(unsigned int), 1, fp);
		unsigned int nsize = Size();
		fwrite(&nsize , sizeof(unsigned int), 1, fp);
		fwrite(m_pd   , sizeof(T), m_nsize, fp);
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
		long			lpos;	// file position of size field
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
	IOResult read(char&   c) { int nr = (int) fread(&c, sizeof(char  ), 1, m_fp); if (nr != 1) return IO_ERROR; return IO_OK; }
	IOResult read(int&    n) { int nr = (int) fread(&n, sizeof(int   ), 1, m_fp); if (nr != 1) return IO_ERROR; return IO_OK; }
	IOResult read(bool&   b) { int nr = (int) fread(&b, sizeof(bool  ), 1, m_fp); if (nr != 1) return IO_ERROR; return IO_OK; }
	IOResult read(float&  f) { int nr = (int) fread(&f, sizeof(float ), 1, m_fp); if (nr != 1) return IO_ERROR; return IO_OK; }
	IOResult read(double& g) { int nr = (int) fread(&g, sizeof(double), 1, m_fp); if (nr != 1) return IO_ERROR; return IO_OK; }

	IOResult read(unsigned int& n) { size_t nr = fread(&n, sizeof(unsigned int), 1, m_fp); if (nr != 1) return IO_ERROR; return IO_OK; }

	IOResult read(char*   pc, int n) { int nr = (int) fread(pc, sizeof(char  ), n, m_fp); if (nr != n) return IO_ERROR; return IO_OK; }
	IOResult read(int*    pi, int n) { int nr = (int) fread(pi, sizeof(int   ), n, m_fp); if (nr != n) return IO_ERROR; return IO_OK; }
	IOResult read(bool*   pb, int n) { int nr = (int) fread(pb, sizeof(bool  ), n, m_fp); if (nr != n) return IO_ERROR; return IO_OK; }
	IOResult read(float*  pf, int n) { int nr = (int) fread(pf, sizeof(float ), n, m_fp); if (nr != n) return IO_ERROR; return IO_OK; }
	IOResult read(double* pg, int n) { int nr = (int) fread(pg, sizeof(double), n, m_fp); if (nr != n) return IO_ERROR; return IO_OK; }

	IOResult read(char* sz)
	{
		IOResult ret;
		int l, nr;
		ret = read(l); if (ret != IO_OK) return ret;
		nr = (int) fread(sz, 1, l, m_fp); if (nr != l) return IO_ERROR;
		sz[l] = 0;
		return IO_OK;
	}

protected:
	FILE*	m_fp;		// the file pointer
	bool	m_bSaving;	// read or write mode?

	// write data
	OBranch*	m_pRoot;	// chunk tree root
	OBranch*	m_pChunk;	// current chunk

	// read data
	bool			m_bend;		// chunk end flag
	stack<CHUNK*>	m_Chunk;
};
