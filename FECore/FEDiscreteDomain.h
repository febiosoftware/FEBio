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
#include "FEDomain.h"

//-----------------------------------------------------------------------------
//! domain for discrete elements
class FECORE_API FEDiscreteDomain : public FEDomain
{
	FECORE_SUPER_CLASS(FEDISCRETEDOMAIN_ID)
	FECORE_BASE_CLASS(FEDiscreteDomain)

public:
	FEDiscreteDomain(FEModel* fem) : FEDomain(FE_DOMAIN_DISCRETE, fem) {}

	bool Create(int nsize, FE_Element_Spec espec) override;
	int Elements() const override { return (int) m_Elem.size(); }
	FEElement& ElementRef(int n) override { return m_Elem[n]; }
	const FEElement& ElementRef(int n) const override { return m_Elem[n]; }

	FEDiscreteElement& Element(int n) { return m_Elem[n]; }

	bool Init() override;

	void Reset() override;

	//! copy data from another domain (overridden from FEDomain)
	void CopyFrom(FEMeshPartition* pd) override;

public:
	void AddElement(int eid, int n[2]);

protected:
	vector<FEDiscreteElement>	m_Elem;
};
