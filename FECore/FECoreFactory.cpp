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



#include "stdafx.h"
#include "FECoreFactory.h"
#include "FECoreBase.h"
#include "log.h"
#include <assert.h>

//-----------------------------------------------------------------------------
//! constructor
FECoreFactory::FECoreFactory(SUPER_CLASS_ID scid, const char* szclass, const char* szbase, const char* szalias, int nspec) : m_scid(scid)
{ 
	m_szclass = szclass;
	m_szalias = szalias;
	m_szbase = szbase;
	m_module = 0; 
	m_spec = nspec;
	m_alloc_id = -1;
}

//-----------------------------------------------------------------------------
//! virtual constructor
FECoreFactory::~FECoreFactory(){}

//-----------------------------------------------------------------------------
//! set the module ID
void FECoreFactory::SetModuleID(unsigned int mid)
{
	m_module = mid;
}

//-----------------------------------------------------------------------------
FECoreBase* FECoreFactory::CreateInstance(FEModel* pfem) const
{
	// create a new instance of this class
	FECoreBase* pclass = Create(pfem); assert(pclass);
	if (pclass == 0) return 0;

	// store the factory that created the class
	pclass->SetFactoryClass(this);

	// build the class descriptor
	if (pclass->BuildClass() == false) return 0;

	if (m_spec != -1)
	{
		int n1 = FECORE_SPEC_MAJOR(m_spec);
		int n2 = FECORE_SPEC_MINOR(m_spec);
		if (pfem) feLogWarningEx(pfem, "\"%s\" is deprecated in spec %d.%d!", m_szalias, n1, n2);
	}

	// return the pointer
	return pclass;
}
