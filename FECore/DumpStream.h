/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in 
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



#pragma once
#include <vector>
#include <string>
#include <map>
#include <string.h>
#include "vec3d.h"
#include "mat3d.h"
#include "quatd.h"
#include "tens3d.h"
#include "fecore_api.h"
#include "FECoreKernel.h"

//-----------------------------------------------------------------------------
class FEModel;

//-----------------------------------------------------------------------------
//! A dump stream is used to serialize data to and from a data stream.
//! This is used in FEBio for running and cold restarts. 
//! This is just an abstract base class. Classes must be derived from this
//! to implement the actual storage mechanism.
class FECORE_API DumpStream
{
	struct Pointer {
		void*			pd;
		unsigned int	id;
	};

public:
	// This class is thrown when an error occurs reading the dumpfile
	class ReadError{};

public:
	//! constructor
	DumpStream(FEModel& fem);

	//! destructor
	virtual ~DumpStream();

	//! See if the stream is used for input or output
	bool IsSaving() const;

	//! See if the stream is used for input
	bool IsLoading() const;

	//! See if shallow flag is set
	bool IsShallow() const;

	// open the stream
	virtual void Open(bool bsave, bool bshallow);

	// get the FE model
	FEModel& GetFEModel() { return m_fem; }

public:
	// These functions must be overloaded by derived classes
	// and should return the number of bytes serialized
	virtual size_t write(const void* pd, size_t size, size_t count) = 0;
	virtual size_t read(void* pd, size_t size, size_t count) = 0;

	// additional function that need to be overridden
	virtual void clear() = 0;

	void check();

	void LockPointerTable();
	void UnlockPointerTable();

public:
	// input-output operators (will call correct operator depending on input or output mode)
	template <typename T> DumpStream& operator & (T& o);

public: // output operators
	DumpStream& operator << (const char* sz);
	DumpStream& operator << (char* sz);
	DumpStream& operator << (const double a[3][3]);
	DumpStream& operator << (std::string& s);
	DumpStream& operator << (const std::string& s);
	DumpStream& operator << (bool b);
	DumpStream& operator << (int n);
	template <typename T> DumpStream& operator << (T& o);
	template <typename T> DumpStream& operator << (std::vector<T>& o);
	template <typename A, typename B> DumpStream& operator << (std::map<A, B>& o);
	template <typename T, std::size_t N> DumpStream& operator << (T(&a)[N]);
	template <typename T> DumpStream& operator << (T* &a);

	template <typename T> DumpStream& write_raw(const T& o);

public: // input operators
	DumpStream& operator >> (char* sz);
	DumpStream& operator >> (double a[3][3]);
	DumpStream& operator >> (std::string& s);
	DumpStream& operator >> (bool& b);
	template <typename T> DumpStream& operator >> (T& o);
	template <typename T> DumpStream& operator >> (std::vector<T>& o);
	template <typename A, typename B> DumpStream& operator >> (std::map<A, B>& o);
	template <typename T, std::size_t N> DumpStream& operator >> (T(&a)[N]);
	template <typename T> DumpStream& operator >> (T* &a);

	template <typename T> DumpStream& read_raw(T& o);

private:
	int FindPointer(void* p);
	int FindPointer(int id);
	void AddPointer(void* p);

private:
	bool		m_bsave;	//!< true if output stream, false for input stream
	bool		m_bshallow;	//!< if true only shallow data needs to be serialized
	FEModel&	m_fem;		//!< the FE Model that is being serialized

	size_t	m_bytes_serialized;	//!< number or bytes serialized

	bool					m_ptr_lock;
	std::vector<Pointer>	m_ptr;
};

template <typename T> DumpStream& DumpStream::write_raw(const T& o)
{
	m_bytes_serialized += write(&o, sizeof(T), 1);
	return *this;
}

template <typename T> DumpStream& DumpStream::read_raw(T& o)
{
	m_bytes_serialized += read(&o, sizeof(T), 1);
	return *this;
}

template <typename T> inline DumpStream& DumpStream::operator & (T& o)
{
	if (IsSaving()) (*this) << o; else (*this) >> o;
	return *this;
}

template <typename T> inline DumpStream& DumpStream::operator << (T& o)
{
	AddPointer((void*)&o);
	o.Serialize(*this);
	check();
	return *this;
}

template <> inline DumpStream& DumpStream::operator << (int&          o) { return write_raw(o); }
template <> inline DumpStream& DumpStream::operator << (unsigned int& o) { return write_raw(o); }
template <> inline DumpStream& DumpStream::operator << (double&   o) { return write_raw(o); }
template <> inline DumpStream& DumpStream::operator << (vec2d&    o) { return write_raw(o); }
template <> inline DumpStream& DumpStream::operator << (vec3d&    o) { return write_raw(o); }
template <> inline DumpStream& DumpStream::operator << (quatd&    o) { return write_raw(o); }
template <> inline DumpStream& DumpStream::operator << (mat2d&    o) { return write_raw(o); }
template <> inline DumpStream& DumpStream::operator << (mat3d&    o) { return write_raw(o); }
template <> inline DumpStream& DumpStream::operator << (mat3ds&   o) { return write_raw(o); }
template <> inline DumpStream& DumpStream::operator << (mat3dd&   o) { return write_raw(o); }
template <> inline DumpStream& DumpStream::operator << (mat3da&   o) { return write_raw(o); }
template <> inline DumpStream& DumpStream::operator << (tens3ds&  o) { return write_raw(o); }
template <> inline DumpStream& DumpStream::operator << (tens3drs& o) { return write_raw(o); }

template <typename T> inline DumpStream& DumpStream::operator >> (T& o)
{
	AddPointer((void*)&o);
	o.Serialize(*this);
	check();
	return *this;
}

template <> inline DumpStream& DumpStream::operator >> (int&          o) { return read_raw(o); }
template <> inline DumpStream& DumpStream::operator >> (unsigned int& o) { return read_raw(o); }
template <> inline DumpStream& DumpStream::operator >> (double&   o) { return read_raw(o); }
template <> inline DumpStream& DumpStream::operator >> (vec2d&    o) { return read_raw(o); }
template <> inline DumpStream& DumpStream::operator >> (vec3d&    o) { return read_raw(o); }
template <> inline DumpStream& DumpStream::operator >> (quatd&    o) { return read_raw(o); }
template <> inline DumpStream& DumpStream::operator >> (mat2d&    o) { return read_raw(o); }
template <> inline DumpStream& DumpStream::operator >> (mat3d&    o) { return read_raw(o); }
template <> inline DumpStream& DumpStream::operator >> (mat3ds&   o) { return read_raw(o); }
template <> inline DumpStream& DumpStream::operator >> (mat3dd&   o) { return read_raw(o); }
template <> inline DumpStream& DumpStream::operator >> (mat3da&   o) { return read_raw(o); }
template <> inline DumpStream& DumpStream::operator >> (tens3ds&  o) { return read_raw(o); }
template <> inline DumpStream& DumpStream::operator >> (tens3drs& o) { return read_raw(o); }

template <typename T> inline DumpStream& DumpStream::operator << (std::vector<T>& o)
{
	int N = (int) o.size();
	write(&N, sizeof(int), 1);
	for (int i=0; i<N; ++i) (*this) << o[i];
	return *this;
}

template <typename T> inline DumpStream& DumpStream::operator >> (std::vector<T>& o)
{
	DumpStream& This = *this;
	int N;
	read(&N, sizeof(int), 1);
	if (N > 0)
	{
		o.resize(N);
		for (int i = 0; i<N; ++i) (*this) >> o[i];
	}
	return This;
}

template <> inline DumpStream& DumpStream::operator << (std::vector<bool>& o)
{
	DumpStream& This = *this;
	int N = (int) o.size();
	write(&N, sizeof(int), 1);
	for (int i=0; i<N; ++i) 
	{
		bool b = o[i];
		This << b;
	}
	return This;
}

template <> inline DumpStream& DumpStream::operator >> (std::vector<bool>& o)
{
	DumpStream& This = *this;
	int N;
	read(&N, sizeof(int), 1);
	if (N > 0)
	{
		o.resize(N);
		for (int i=0; i<N; ++i) 
		{
			bool b;
			This >> b;
			o[i] = b;
		}
	}
	return This;
}

template <typename A, typename B> DumpStream& DumpStream::operator << (std::map<A, B>& o)
{
	DumpStream& ar = *this;
	int N = (int)o.size();
	ar << N;
    for (typename std::map<A, B>::iterator it = o.begin(); it != o.end(); ++it)
	{
		const A& a = it->first;
		B& b = it->second;
		ar << a << b;
	}
	return ar;
}

template <typename A, typename B> DumpStream& DumpStream::operator >> (std::map<A, B>& o)
{
	DumpStream& ar = *this;
	int N = 0;
	ar >> N;
	o.clear();
	for (int i=0; i<N; ++i)
	{
		A a;
		B b;
		ar >> a >> b;
		o[a] = b;
	}
	return ar;
}

template <typename T, std::size_t N> DumpStream& DumpStream::operator << (T(&a)[N])
{
	for (int i = 0; i < N; ++i) (*this) << a[i];
	return *this;
}

template <typename T, std::size_t N> DumpStream& DumpStream::operator >> (T(&a)[N])
{
	for (int i = 0; i < N; ++i) (*this) >> a[i];
	return *this;
}

template <typename T> DumpStream& DumpStream::operator << (T* &a)
{
	DumpStream& ar = *this;

	// see if we already stored this pointer 
	int pid = FindPointer((void*)a);
	ar << pid;
	if (pid != -1) return ar;

	// store the pointer in the table
	AddPointer((void*)a);

	// If we are storing a deep copy we need to store class info 
	// so that we can reinstantiate the class
	if (ar.IsShallow() == false)
	{
		// store the class info
		T::SaveClass(*this, a);
	}

	// serialize the object (assuming it has a Serialize member)
	a->Serialize(*this);

	return *this;
}

template <typename T> DumpStream& DumpStream::operator >> (T* &a)
{
	DumpStream& ar = *this;

	// get the pointer id
	int pid;
	ar >> pid;
	if (pid != -1)
	{
		a = (T*)(m_ptr[pid].pd);
		return ar;
	}

	// read class identifier and instatiate class
	if (ar.IsShallow() == false)
	{
		a = dynamic_cast<T*>(T::LoadClass(ar, a));
	}

	// store the pointer
	AddPointer((void*)a);

	// serialize the object
	a->Serialize(*this);

	return *this;
}
