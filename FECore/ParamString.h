/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
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
#include "fecore_api.h"

//-----------------------------------------------------------------------------
class ParamRef
{
public:
	int			_index;		// zero-based index of parameter (-1 if not available)
	int			_id;		// ID of parameter (-1 if not available)
	std::string	_name;		// name of parameter
	std::string	_idName;	// string ID of parameter

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
class FECORE_API ParamString
{
public:
	//! constructor
	ParamString(const char* sz);

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
	bool operator == (const std::string& s) const;

	//! compare to a string
	bool operator != (const std::string& s) const;

	//! Get the ID (-1 if ID not a number)
	int ID() const;

	//! get the index (-1 if index not used)
	int Index() const;

	//! get the index name (null if not defined)
	const char* IDString() const;

	//! get the zero-valued string 
	const char* c_str() const;

	//! return a string
	std::string string() const;

private:
	ParamString() {}

private:
	std::vector<ParamRef>	m_param;
};
