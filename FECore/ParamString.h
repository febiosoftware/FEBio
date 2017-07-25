#pragma once
#include <vector>
#include <string>
#include "fecore_export.h"
using namespace std;

//-----------------------------------------------------------------------------
class FECOREDLL_EXPORT ParamRef
{
public:
	string	_name;		// name of parameter
	int		_index;		// zero-based index of parameter (-1 if not available)
	int		_id;		// ID of parameter (-1 if not available)
	string	_idName;	// string ID of parameter

public:
	ParamRef() : _id(-1), _index(-1) {}
	ParamRef(const ParamRef& p)
	{
		_name = p._name;
		_id = p._id;
		_index = p._index;
		_idName = p._idName;
	}

	void operator = (const ParamRef& p)
	{
		_name = p._name;
		_id = p._id;
		_index = p._index;
		_idName = p._idName;
	}

	void clear()
	{
		_name.clear();
		_id = _index = -1;
		_idName.clear();
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

	//! is the string valid
	bool isValid() const;

	//! clear the string
	void clear();

public:

	//! compare to a string
	bool operator == (const string& s) const;

	//! compare to a string
	bool operator != (const string& s) const;

	//! Get the ID (-1 if ID not a number)
	int ID() const;

	//! get the index (-1 if index not used)
	int Index() const;

	//! get the index name (null if not defined)
	const char* IDString() const;

	//! get the zero-valued string 
	const char* c_str() const;

private:
	ParamString() {}

private:
	vector<ParamRef>	m_param;
};
