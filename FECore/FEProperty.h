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
#include "fecore_api.h"
#include "fecore_enum.h"

//-----------------------------------------------------------------------------
class FECoreBase;
class DumpStream;

//-----------------------------------------------------------------------------
//! A property of a class reflects a member variable of the class that is a 
//! pointer to a FECoreBase derived class. 
class FECORE_API FEProperty
{
public:
	enum Flags
	{
		Optional		= 0x00,
		Required		= 0x01,		// the property is required (default)
		Preferred		= 0x02,		// the property is not required, but a default should be allocated when possible.
		Reference       = 0x04,		// references another class in the model
		Fixed			= 0x08,		// fixed properties are fixed type class members
		TopLevel		= 0x10,		// This is a "top-level" property. 
	};

private:
	//! Name of the property.
	//! Note that the name is not copied so it must point to a static string.
	const char*		m_szname;
	const char*		m_szlongname;	// long name (optional; used in FEBio Studio)
	unsigned int	m_flags;		// bitwise or of flags defined above 

	const char* m_szdefaultType;	// default type string (used by FEBio Studio to initialize required properties).

protected:
	const char* m_className;	// name of class that can be assigned to this

public:
	// Set\Get the name of the property
	FEProperty& SetName(const char* sz);
	const char* GetName() const;

	// get\set the long name
	FEProperty& SetLongName(const char* sz);
	const char* GetLongName() const;

	// get the class name
	const char* GetClassName() const { return m_className; }

	// is the property required
	bool IsRequired() const { return (m_flags & Required) != 0; }

	// is the property preferred
	bool IsPreferred() const { return (m_flags & Preferred) != 0; }

	// is this a reference property
	bool IsReference() const { return (m_flags & Reference) != 0; }

	// is this a top-level property
	bool IsTopLevel() const { return (m_flags & TopLevel) != 0; }

	// set the flags
	void SetFlags(unsigned int flags) { m_flags = flags; }

	// add a flag
	void AddFlag(unsigned int flag) { m_flags |= flag; }

	// get the flags
	unsigned int Flags() const { return m_flags; }

	// get default type (can be null)
	const char* GetDefaultType() const;

	// set the default type
	FEProperty& SetDefaultType(const char* szdefType);

public: // these functions have to be implemented by derived classes

	//! helper function for identifying if this is an array property or not
	virtual bool IsArray() const = 0;

	//! see if the pc parameter is of the correct type for this property
	virtual bool IsType(FECoreBase* pc) const = 0;

	//! set the property
	virtual void SetProperty(FECoreBase* pc) = 0;

	//! return the size of the property
	virtual int size() const = 0;

	//! return a specific property by index
	virtual FECoreBase* get(int i) = 0;

	//! return a specific property by name
	virtual FECoreBase* get(const char* szname) = 0;

	//! return a specific property by ID
	virtual FECoreBase* getFromID(int nid) = 0;

	//! serialize property data
	virtual void Serialize(DumpStream& ar) = 0;

	//! initializatoin
	virtual bool Init() = 0;

	//! validation
	virtual bool Validate() = 0;

	//! Get the parent of this property
	FECoreBase* GetParent() { return m_pParent; }

	//! Set the parent of this property
	virtual void SetParent(FECoreBase* parent) { m_pParent = parent; }

	//! Get the class ID
	SUPER_CLASS_ID GetSuperClassID() const { return m_superClassID; }

public:
	virtual ~FEProperty();

protected:
	//! some helper functions for reading, writing properties
	void Write(DumpStream& ar, FECoreBase* pc);
	FECoreBase* Read(DumpStream& ar);

protected:
	// This class should not be created directly
	FEProperty(SUPER_CLASS_ID classID);

protected:
	FECoreBase*		m_pParent;		//!< pointer to the parent class (i.e. the class that defines this property)
	SUPER_CLASS_ID	m_superClassID;	//!< The super class ID
};
