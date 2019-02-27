#pragma once
#include <vector>
#include <string>
#include <string.h>
#include "vec3d.h"
#include "mat3d.h"
#include "quatd.h"
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
	virtual size_t write(const void* pd, size_t size, size_t count) = 0;
	virtual size_t read(void* pd, size_t size, size_t count) = 0;
	virtual void clear() = 0;
	virtual void check(){}

public:
	// input-output operators (will call correct operator depending on input or output mode)
	template <typename T> DumpStream& operator & (T& o);

public: // output operators
	DumpStream& operator << (const char* sz);
	DumpStream& operator << (char* sz);
	DumpStream& operator << (const double a[3][3]);
	DumpStream& operator << (const std::string& s);
	template <typename T> DumpStream& operator << (const T& o);
	template <typename T> DumpStream& operator << (std::vector<T>& o);
	template <typename T, std::size_t N> DumpStream& operator << (T(&a)[N]);
	template <typename T> DumpStream& operator << (T* const &a);

public: // input operators
	DumpStream& operator >> (char* sz);
	DumpStream& operator >> (double a[3][3]);
	DumpStream& operator >> (std::string& s);
	template <typename T> DumpStream& operator >> (T& o);
	template <typename T> DumpStream& operator >> (std::vector<T>& o);
	template <typename T, std::size_t N> DumpStream& operator >> (T(&a)[N]);
	template <typename T> DumpStream& operator >> (T* &a);

private:
	bool		m_bsave;	//!< true if output stream, false for input stream
	bool		m_bshallow;	//!< if true only shallow data needs to be serialized
	FEModel&	m_fem;		//!< the FE Model that is being serialized
};

template <typename T> inline DumpStream& DumpStream::operator & (T& o)
{
	if (IsSaving()) (*this) << o; else (*this) >> o;
	return *this;
}

template <typename T> inline DumpStream& DumpStream::operator << (const T& o)
{
	write(&o, sizeof(T), 1);
	return *this;
}

template <typename T> inline DumpStream& DumpStream::operator >> (T& o)
{
	read(&o, sizeof(T), 1);
	return *this;
}

template <> inline DumpStream& DumpStream::operator << (const bool& o)
{
	int b = (o==true?1:0);
	write(&b, sizeof(int), 1);
	return *this;
}

template <> inline DumpStream& DumpStream::operator >> (bool& o)
{
	int b;
	read(&b, sizeof(int), 1);
	o = (b==1);
	return *this;
}

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
		for (int i=0; i<N; ++i) 
		{
			bool b;
			This >> b;
			o[i] = b;
		}
	}
	return This;
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

template <typename T> DumpStream& DumpStream::operator << (T* const &a)
{
	// store the class identifier
	int classID = 0;
	if (a == nullptr) { (*this) << classID; return *this; }
	classID = a->GetSuperClassID();
	const char* sztype = a->GetTypeStr();
	(*this) << classID;
	(*this) << sztype;

	// serialize the object (assuming it has a Serialize member)
	a->Serialize(*this);

	return *this;
}

template <typename T> DumpStream& DumpStream::operator >> (T* &a)
{
	DumpStream& ar = *this;
	// read class identifier and instatiate class
	int classID = 0;
	ar >> classID;
	if (classID == FEINVALID_ID) {
		a = nullptr; return ar;
	}
	char sztype[256] = { 0 };
	ar >> sztype;

	// instantiate the class
	a = fecore_new<T>(classID, sztype, &GetFEModel());
	assert(a);
	if (a == nullptr) throw ReadError();

	// serialize the object
	a->Serialize(*this);

	return *this;
}
