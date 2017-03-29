#pragma once
#include <vector>
#include <string>
#include "fecore_export.h"
using namespace std;

//-----------------------------------------------------------------------------
class FECOREDLL_EXPORT ParamRef
{
public:
	string	name;		// name of parameter
	int		id;			// index of parameter (-1 if not available)
	string	idName;		// string index

public:
	ParamRef() : id(-1){}
	ParamRef(const ParamRef& p)
	{
		name = p.name;
		id = p.id;
		idName = p.idName;
	}

	void operator = (const ParamRef& p)
	{
		name = p.name;
		id = p.id;
		idName = p.idName;
	}
};

//-----------------------------------------------------------------------------
// helper class for retrieving parameters
class FECOREDLL_EXPORT ParamString
{
public:
	//! constructor
	explicit ParamString(const char* sz);

	//! copy constructor
	ParamString(const ParamString& p);

	//! assignment operator
	void operator = (const ParamString& p);

	//! destructor
	~ParamString();

	//! number of refs in string
	int count() const;

	//! return a new string starting from the next component
	ParamString next() const;

	//! return the last string
	ParamString last() const;

public:

	//! compare to a string
	bool operator == (const string& s) const;

	//! compare to a string
	bool operator != (const string& s) const;

	//! get the index (-1 if index not a number)
	int Index() const;

	//! get the index name (null if not defined)
	const char* IndexString() const;

	//! get the zero-valued string 
	const char* c_str() const;

private:
	ParamString() {}

private:
	vector<ParamRef>	m_param;
};
