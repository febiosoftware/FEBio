// Archive.h: interface for the Archive class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_ARCHIVE_H__B95A81B1_BBFB_46E5_B9B3_7675ED8A6029__INCLUDED_)
#define AFX_ARCHIVE_H__B95A81B1_BBFB_46E5_B9B3_7675ED8A6029__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <stdio.h>
#include <string.h>
#include "vec3d.h"
#include "mat3d.h"
#include "quatd.h"
#include "FECore/vector.h"

//-----------------------------------------------------------------------------
//! Class for serializing data to a binary archive.

//! This class is used to read data from or write
//! data to a binary file. The class defines several operators to 
//! simplify in- and output.
//! \sa FEM::Serialize()

class Archive  
{
public:
	Archive();
	virtual ~Archive();

	//! Open archive for reading
	bool Open(const char* szfile);

	//! Open archive for writing
	bool Create(const char* szfile);

	//! Open archive for appending
	bool Append(const char* szfile);

	//! Close archive
	void Close();

	//! Check mode
	bool IsSaving() { return m_bsave; }

	//! See if the archive is valid
	bool IsValid() { return (m_fp != 0); }

	//! Flush the archive
	void Flush() { fflush(m_fp); }

	//@{ 
	//! output operators
	Archive& operator << (const char* sz) 
	{ 
		int n = strlen(sz); 
		fwrite(&n, sizeof(int), 1, m_fp);
		fwrite(sz, sizeof(char), n, m_fp);
		return (*this);
	}
	Archive& operator << (const double a[3][3])
	{
		fwrite(a, sizeof(double), 9, m_fp);
		return (*this);
	}

	template <class T> Archive& operator << (const T& o) { fwrite(&o, sizeof(T), 1, m_fp); return (*this); }

	template <class T> Archive& operator << (vector<T>& v)
	{
		int n = v.size();
		fwrite(&n, sizeof(int), 1, m_fp);
		fwrite((T*) &v[0], sizeof(T), v.size(), m_fp);
		return (*this);
	}
	//@}


	//@{
	//! input operators
	Archive& operator >> (char* sz) 
	{ 
		int n;
		fread(&n, sizeof(int), 1, m_fp);
		fread(sz, sizeof(char), n, m_fp);
		sz[n] = 0;
		return (*this);
	}

	Archive& operator >> (double a[3][3])
	{
		fread(a, sizeof(double), 9, m_fp);
		return (*this);
	}

	template <class T> Archive& operator >> (T& o) { fread(&o, sizeof(T), 1, m_fp); return (*this); }

	template <class T> Archive& operator >> (vector<T>& v)
	{
		int n;
		fread(&n, sizeof(int), 1, m_fp);
		v.resize(n);
		fread((T*) &v[0], sizeof(T), n, m_fp);
		return (*this);
	}
	//@}


	//! write buffer to archive
	size_t write(const void* pd, size_t size, size_t count)
	{
		return fwrite(pd, size, count, m_fp);
	}

	//! read buffer from archive
	size_t read(void* pd, size_t size, size_t count)
	{
		return fread(pd, size, count, m_fp);
	}


protected:
	FILE*	m_fp;		//!< The actual file pointer
	bool	m_bsave;	//!< Save flag
};

#endif // !defined(AFX_ARCHIVE_H__B95A81B1_BBFB_46E5_B9B3_7675ED8A6029__INCLUDED_)
