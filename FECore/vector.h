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
#include <vector>
#include <algorithm>
using namespace std;

double operator*(const vector<double>& a, const vector<double>& b);
template<typename T> void zero(vector<T>& a) { fill(a.begin(), a.end(), T(0)); }
template<typename T> void set(vector<T>& a, const T& v) { fill(a.begin(), a.end(), v); }
vector<double>& operator += (vector<double>& a, const vector<double>& b);

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

	void create(int n) 
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
	~ptr_vector() { clear(); }

	// operators
	T& operator [] (int i) { return *m_pdata[i]; }
	const T& operator [] (int i) const { return *m_pdata[i]; }

	// function to return the size of the vector
	int size() const { return m_nsize; }

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
