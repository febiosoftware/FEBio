// Archive.h: interface for the Archive class.
//
//////////////////////////////////////////////////////////////////////
#pragma once
#include <vector>
#include <queue>
#include <string.h>
using namespace std;

//-----------------------------------------------------------------------------
//! Class used for writing data to a binary file
class Archive
{
public:
	// con/de-structor
	Archive();
	~Archive();

	//! close the file
	void Close();

	//! create the file
	bool Create(const char* sz);

	//! create a branch
	void BeginChunk(unsigned int id);

	//! end branch
	void EndChunk();

	//! write a single data entry
	template <typename T> void WriteChunk(unsigned int id, const T& d)
	{
		unsigned int n[3] = {id, sizeof(T), 1};
		fwrite(n, sizeof(unsigned int), 3, m_fp);
		fwrite(&d, sizeof(T), 1, m_fp);
	}

	//! write a vector of data
	template <typename T> void WriteChunk(unsigned int id, const vector<T>& v)
	{
		unsigned int nv = (unsigned int) v.size();
		unsigned int n[3] = {id, sizeof(T), nv};
		fwrite(n, sizeof(unsigned int), 3, m_fp);
		if (nv > 0) fwrite(&v[0], sizeof(T), nv, m_fp);
	}

	//! write array data
	template <typename T> void WriteChunk(unsigned int id, const T* pv, size_t nv)
	{
		unsigned int n[3] = {id, sizeof(T), nv};
		fwrite(n, sizeof(unsigned int), 3, m_fp);
		if (nv > 0) fwrite(pv, sizeof(T), nv, m_fp);
	}

	void WriteChunk(unsigned int id, const char* sz)
	{
		int l = strlen(sz);
		unsigned int n[3] = {id, sizeof(char), l};
		if (l > 0) fwrite(sz, sizeof(char), l, m_fp);
	}

	template <typename T> Archive& operator << (const T& d)
	{
		fwrite(&d, sizeof(T), 1, m_fp);
		return (*this);
	}

private:
	FILE*	m_fp;
	queue<unsigned int>	m_tag;
};
