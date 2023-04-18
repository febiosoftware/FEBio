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
#include "FEDomainList.h"
#include "DumpStream.h"
#include "FEDomain.h"
#include "FEMesh.h"
#include <assert.h>
using namespace std;

FEDomainList::FEDomainList()
{

}

FEDomainList::FEDomainList(FEDomainList& domList)
{
	m_dom = domList.m_dom;
}

//! Clear the domain list
void FEDomainList::Clear()
{
	m_dom.clear();
}

void FEDomainList::AddDomain(FEDomain* dom)
{
	// see if this domain is already a member of this list
	if (IsMember(dom))
	{
//		assert(false);
		return;
	}

	// it's not, so let's add it
	m_dom.push_back(dom);
}

//! Add a domain list
void FEDomainList::AddDomainList(const FEDomainList& domList)
{
	for (int i = 0; i < domList.Domains(); ++i)
	{
		FEDomain* d = const_cast<FEDomain*>(domList.GetDomain(i));
		AddDomain(d);
	}
}

bool FEDomainList::IsMember(const FEDomain* dom) const
{
	// loop over all the domains
	for (size_t i = 0; i < m_dom.size(); ++i)
	{
		if (m_dom[i] == dom)
		{
			// found it!
			return true;
		}
	}

	// better luck next time!
	return false;
}

//! serialization
void FEDomainList::Serialize(DumpStream& ar)
{
	if (ar.IsShallow()) return;
	ar & m_dom;
}
