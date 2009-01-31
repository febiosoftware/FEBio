// vector.h: interface for the vector class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_VECTOR_H__9F132D73_20B9_4AE9_A40B_EE4FB9D0FABD__INCLUDED_)
#define AFX_VECTOR_H__9F132D73_20B9_4AE9_A40B_EE4FB9D0FABD__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <math.h>
#include <memory.h>

//-----------------------------------------------------------------------------
//! This class implements a dynamic array.

//! We do not use the STL vector class because it does not allow us to create
//! vectors of arrays, e.g. vector<int[3]>. Our class does allow such a cool
//! construct!

template <typename T> class vector  
{
public:
	// constructors
	vector() : m_pdata(0), m_nsize(0), m_nbufsize(0), m_ngrowsize(1) {}
	explicit vector(int n) : m_pdata(new T[n]), m_nsize(n), m_nbufsize(n), m_ngrowsize(1) {}
	explicit vector(int n, int nalloc) : m_ngrowsize(1)
	{
		m_nsize = 0;
		m_nbufsize = 0;
		m_pdata = 0;
		if (nalloc < n) nalloc = n;
		alloc(nalloc);
		m_nsize = n; 
	}
	vector(const vector<T>& a);

	virtual void create(int n, int nalloc=0) 
	{ 
		if (nalloc < n) nalloc = n;
		alloc(nalloc);
		m_nsize = n; 
	}

	// destructor
	virtual ~vector() { delete [] m_pdata; m_pdata = 0; }

	void clear()
	{
		m_nsize = 0;
	}

	// operators
	T& operator [] (int i) { return m_pdata[i]; }
	const T& operator [] (int i) const { return m_pdata[i]; }

	vector<T> operator * (double g)
	{
		vector<T> r(m_nsize);
		for (int i=0; i<m_nsize; ++i) r.m_pdata[i] = m_pdata[i]*g;
		return r;
	}

	vector<T>& operator *= (const double g)
	{
		for (int i=0; i<m_nsize; ++i) m_pdata[i] *= g;
		return (*this);
	}

	vector<T>& operator -= (vector<T>& a)
	{
		for (int i=0; i<m_nsize; ++i) m_pdata[i] -= a[i];
		return (*this);
	}

	vector<T>& operator += (const vector<T>& a)
	{
		for (int i=0; i<m_nsize; ++i) m_pdata[i] += a[i];
		return (*this);
	}

	vector<T>& operator = (const vector<T>& a);
	operator T*
		() { return m_pdata; }

	vector<T> operator - (const vector<T>& a) const
	{
		vector<T> b(m_nsize);
		for (int i=0; i<m_nsize; ++i) b.m_pdata[i] = m_pdata[i] - a.m_pdata[i];

		return b;
	}

	vector<T> operator + (const vector<T>& a) const
	{
		vector<T> b(m_nsize);
		for (int i=0; i<m_nsize; ++i) b.m_pdata[i] = m_pdata[i] + a.m_pdata[i];

		return b;
	}

	// function to set all vector items to zero
	// note that the conversion '(T) 0' must make sense when using this function
	void zero() { memset(m_pdata, 0, m_nsize*sizeof(T)); }

	// function to set all items to a certain value
	void set(const T& a) { for (int i=0; i<m_nsize; ++i) m_pdata[i] = a; }

	// function to return the size of the vector
	int size() const { return m_nsize; }

	// append operation
	// use only when copy operation is valid for class T
	void add(T& a)
	{
		if (m_nsize + 1 > m_nbufsize) realloc(m_nbufsize + m_ngrowsize);
		
		m_pdata[m_nsize] = a;
		++m_nsize;
	}

	// insert operator
	void insert(T& a, int n)
	{
		T* ptmp = m_pdata;
		if (m_nsize + 1 > m_nbufsize)
		{
			m_nbufsize += m_ngrowsize;
			ptmp = new T[m_nbufsize];
			for (int i=0; i<n; ++i) ptmp[i] = m_pdata[i];
		}

		for (int i=m_nsize-1; i>n; --i) ptmp[i+1] = m_pdata[i];

		ptmp[n] = a;

		if (ptmp != m_pdata)
		{
			delete [] m_pdata;
			m_pdata = ptmp;
		}

		++m_nsize;
	}

	// set the growsize of the vector
	void setgrowsize(int n) { m_ngrowsize = n; }

	// resize array to a specific size
	// use only when copy operation is valid for class T
	void setsize(int n) { realloc(n); m_nsize = n; }

protected:
	// this function allocates a buffer for storage of a specific size
	void alloc(int n)
	{
		if (n > m_nbufsize)
		{
			m_nbufsize = n;
			if (m_pdata) delete [] m_pdata;
			m_pdata = new T[n];
		}
	}

	// this function allocates a buffer and copies (part of) the old data into the new buffer
	void realloc(int n)
	{
		if (n > m_nbufsize)
		{
			T* ptmp = new T[n];
			for (int i=0; i<m_nsize; ++i) ptmp[i] = m_pdata[i];
			if (m_pdata) delete [] m_pdata;
			m_pdata = ptmp;
			m_nbufsize = n;
		}
	}

protected:
	int	m_nsize;		// size of vector
	int	m_nbufsize;		// actual size of allocated buffer
	int	m_ngrowsize;	// size to increase bufsize
	T*	m_pdata;		// the vector data
};

// copy constructor
template <class T> vector<T>::vector(const vector<T>& a)
{
	m_pdata = new T[a.m_nsize];
	m_nsize = a.m_nsize;
	m_nbufsize = m_nsize;
	m_ngrowsize = a.m_ngrowsize;

	for (int i=0; i<m_nsize; ++i) m_pdata[i] = a.m_pdata[i];
}

// assignment operator
template <class T> vector<T>& vector<T>::operator = (const vector<T>& a)
{
	// allocate storage
	create(a.m_nsize);

	// copy data
	for (int i=0; i<m_nsize; ++i) m_pdata[i] = a.m_pdata[i];

	// return reference
	return (*this);
}

double operator*(const vector<double>& a, const vector<double>& b);

///////////////////////////////////////////////////////////////////////////////
// CLASS: ptr_vector
// The ptr_vector class stores an array of pointers. The pointers will be 
// deleted when the vector goes out of scope.
//

template <typename T> class ptr_vector  
{
public:
	// constructors
	ptr_vector() : m_pdata(0), m_nsize(0) {}
	explicit ptr_vector(int n) : m_pdata(new T*[n]), m_nsize(n) { zero(); }

	virtual void create(int n) 
	{ 
		clear();
		m_pdata = new T*[n];
		m_nsize = n; 
		zero();
	}

	void zero() { if (m_pdata) for (int i=0; i<m_nsize; ++i) m_pdata[i] = 0; }

	void add(T* pt)
	{
		T** ptmp = new T*[m_nsize+1];
		for (int i=0; i<m_nsize; ++i) ptmp[i] = m_pdata[i];
		ptmp[m_nsize] = pt;
		delete [] m_pdata;
		m_pdata = ptmp;
		++m_nsize;
	}

	void add(T& pt)
	{
		T** ptmp = new T*[m_nsize+1];
		for (int i=0; i<m_nsize; ++i) ptmp[i] = m_pdata[i];
		ptmp[m_nsize] = new T(pt);
		delete [] m_pdata;
		m_pdata = ptmp;
		++m_nsize;
	}


	void setitem(int n, T* pi) { m_pdata[n] = pi; }

	// destructor
	virtual ~ptr_vector() { clear(); }

	// operators
	T& operator [] (int i) { return *m_pdata[i]; }
	const T& operator [] (int i) const { return *m_pdata[i]; }

	// function to return the size of the vector
	int size() const { return m_nsize; }

protected:
	void clear()
	{
		for (int i=0; i<m_nsize; ++i) delete m_pdata[i];
		delete [] m_pdata; 
		m_pdata = 0;
		m_nsize = 0;
	}

private:
	// do not allow copy construction or assignment operation
	ptr_vector(const ptr_vector<T>& a){}
	void operator = (const ptr_vector<T>& a) {}

protected:
	T**	m_pdata;		// the vector data
	int	m_nsize;		// size of vector
};

#endif // AFX_VECTOR_H__9F132D73_20B9_4AE9_A40B_EE4FB9D0FABD__INCLUDED_
