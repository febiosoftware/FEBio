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
#include "fecore_enum.h"
#include "FEParameterList.h"

//-----------------------------------------------------------------------------
//! Forward declaration of the FEModel class. All classes that register
//! with the framework take a pointer to FEModel as their constructor parameter.
class FEModel;
class FECoreBase;

//-----------------------------------------------------------------------------
//! The factory class contains the mechanism for instantiating a class.
class FECORE_API FECoreFactory : public FEParamContainer
{
public:
	//! constructor
	FECoreFactory(SUPER_CLASS_ID scid, const char* szclass, const char* szbase, const char* szalias, int nspec = -1);

	//! virtual constructor
	virtual ~FECoreFactory();

	//! This is the function that the kernel will use to intantiate an object
	FECoreBase* CreateInstance(FEModel* pfem) const;

public:
	// return the class name
	const char* GetClassName() const { return m_szclass; }

	// return the base class name
	const char* GetBaseClassName() const { return m_szbase; }

	// return the type string identifier
	const char* GetTypeStr() const { return m_szalias; }

	//! return the super-class ID
	SUPER_CLASS_ID GetSuperClassID() const { return m_scid; }

	//! return the module name
	unsigned int GetModuleID() const { return m_module; }

	//! set the module name
	void SetModuleID(unsigned int nid);

	//! Get the spec number
	int GetSpecID() const { return m_spec; }

	//! Set the allocator ID
	void SetAllocatorID(int alloc) { m_alloc_id = alloc; }

	//! Get the allocator ID
	int GetAllocatorID() const { return m_alloc_id; }
	
public:
	//! derived classes implement this to create an instance of a class
	virtual FECoreBase* Create(FEModel*) const = 0;

private:
	const char*		m_szclass;	//!< class name
	const char*		m_szbase;	//!< base class name
	const char*		m_szalias;	//!< class alias string
	int				m_spec;		//!< The max spec number for which this feature is defined (-1 is don't care)
	unsigned int	m_module;	//!< ID of module this class belongs to
	SUPER_CLASS_ID	m_scid;		//!< the super-class ID
	int				m_alloc_id;	//!< allocator ID
};

//-----------------------------------------------------------------------------
//! Forward declarations of classes used by the domain factory
class FEDomain;
class FEMesh;
class FEMaterial;
class FEModel;

//-----------------------------------------------------------------------------
//! Creation of domains are a little more elaborate and deviate from the usual
//! factory methods.
class FECORE_API FEDomainFactory
{
public:
	FEDomainFactory(){}
	virtual ~FEDomainFactory(){}

	virtual FEDomain* CreateDomain(const FE_Element_Spec& spec, FEMesh* pm, FEMaterial* pmat) = 0;
};

#define FECORE_SPEC(major, minor) ((major << 8) + minor)
#define FECORE_SPEC_MAJOR(n) ((n) >> 8)
#define FECORE_SPEC_MINOR(n) ((n) & 0x0F)

// macro for tagging a feature as experimental (i.e. in development)
#define FECORE_EXPERIMENTAL	int(0xFFFF)
