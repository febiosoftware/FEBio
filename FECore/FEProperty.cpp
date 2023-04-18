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
#include "FEProperty.h"
#include "FECoreBase.h"
#include "DumpStream.h"
#include "FECoreKernel.h"
#include "DumpStream.h"

//-----------------------------------------------------------------------------
FEProperty::FEProperty(SUPER_CLASS_ID superClassID) : m_szname(nullptr), m_className(nullptr), m_szlongname(nullptr), m_szdefaultType(nullptr), m_flags(0), m_superClassID(superClassID) {}

//-----------------------------------------------------------------------------
FEProperty::~FEProperty(){}

//-----------------------------------------------------------------------------
//! Set the name of the property.
//! Note that the name is not copied so it must point to a static string.
FEProperty& FEProperty::SetName(const char* sz)
{
	m_szname = sz;
	return *this;
}

//-----------------------------------------------------------------------------
//! Return the name of this property
const char* FEProperty::GetName() const { return m_szname; }

//-----------------------------------------------------------------------------
FEProperty& FEProperty::SetLongName(const char* sz)
{
	m_szlongname = sz;
	return *this;
}

//-----------------------------------------------------------------------------
const char* FEProperty::GetLongName() const { return m_szlongname; }

//-----------------------------------------------------------------------------
const char* FEProperty::GetDefaultType() const
{
	return m_szdefaultType;
}

//-----------------------------------------------------------------------------
FEProperty& FEProperty::SetDefaultType(const char* szdefType)
{
	m_szdefaultType = szdefType;
	return *this;
}

//-----------------------------------------------------------------------------
void FEProperty::Write(DumpStream& ar, FECoreBase* pc)
{
	int nflag = (pc == 0 ? 0 : 1);
	ar << nflag;
	if (nflag)
	{
		int ntype = (int) pc->GetSuperClassID();
		ar << pc->GetTypeStr();
		ar << ntype;
		pc->Serialize(ar);
	}
}

//-----------------------------------------------------------------------------
FECoreBase* FEProperty::Read(DumpStream& ar)
{
	int nflag = 0;
	FECoreBase* pm = 0;
	ar >> nflag;
	if (nflag)
	{
		char sz[256];
		int ntype = FEINVALID_ID;
		ar >> sz;
		ar >> ntype;
		pm = fecore_new<FECoreBase>(ntype, sz, &ar.GetFEModel());
		pm->SetParent(GetParent());
		pm->Serialize(ar);

		// TODO: Do I really need to do this here?
		//pm->Init();
	}
	return pm;
}
