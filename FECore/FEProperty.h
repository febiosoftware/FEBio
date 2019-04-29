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
#include "fecore_api.h"
#include "fecore_enum.h"

//-----------------------------------------------------------------------------
class FEParam;
class FECoreBase;
class ParamString;
class DumpStream;

//-----------------------------------------------------------------------------
//! First attempt at conceptualizing material properties.
//! A material property is essentially any material that is nested inside
//! another material definition. To faciliate automation of properties, we
//! created an explicit interface for such properties.
//! \todo I'd like to make this available for all FECoreBase classes, not just materials.
class FECORE_API FEProperty
{
public:
	enum Flags
	{
		Optional		= 0x00,
		Required		= 0x01,		// the property is required (default)
		Reference       = 0x02		// references another class in the model
	};

private:
	//! Name of the property.
	//! Note that the name is not copied so it must point to a static string.
	const char*		m_szname;
	unsigned int	m_flags;	// true if this flag is required (false if optional). Used in FEMaterial::Init().

public:
	// Set\Get the name of the property
	FEProperty& SetName(const char* sz);
	const char* GetName() const;

	// is the property required
	bool IsRequired() const { return (m_flags & Required) != 0; }

	// is this a referemce property
	bool IsReference() const { return (m_flags & Reference) != 0; }

	// set the flags
	void SetFlags(unsigned int flags) { m_flags = flags; }

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
	SUPER_CLASS_ID GetClassID() const { return m_classID; }

protected:
	//! some helper functions for reading, writing properties
	void Write(DumpStream& ar, FECoreBase* pc);
	FECoreBase* Read(DumpStream& ar);

protected:
	// This class should not be created directly
	FEProperty(SUPER_CLASS_ID classID);
	virtual ~FEProperty();

protected:
	FECoreBase*		m_pParent;	//!< pointer to the "parent" material
	SUPER_CLASS_ID	m_classID;	//!< The class ID
};
