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
#include "fecore_api.h"
#include "FEItemList.h"
#include "FEElement.h"
#include "FEDomainList.h"
#include "FENodeList.h"
#include <vector>
#include <string>

//-----------------------------------------------------------------------------
class FEMesh;
class DumpStream;

//-----------------------------------------------------------------------------
// This class defines a set of elements
class FECORE_API FEElementSet : public FEItemList
{
public:
	//! constructor
	FEElementSet(FEModel* fem);	// TODO: remove!
	FEElementSet(FEMesh* fem);

	// Create the element set
	void Create(const std::vector<int>& elemList);
	void Create(FEDomain* dom, const std::vector<int>& elemList);

	void CopyFrom(FEElementSet& eset);

	// add another element set
	void Add(const FEElementSet& set);

	// Create the element set from a domain
	void Create(FEDomain* dom);

	// Create the element set from a domain
	void Create(FEDomainList& dom);

	// Return number of elements in the set
	int Elements() const { return (int)m_Elem.size(); }

	int operator [] (int i) const { return m_Elem[i]; }

	// return the local index of an element into the element set
	// returns -1 if the element is not part of element set
	int GetLocalIndex(const FEElement& el) const;
	bool Contains(const FEElement& el) const;

	// Get the element ID list
	const std::vector<int>& GetElementIDList() const { return m_Elem; }

	// get the domain list that generated the element set
	FEDomainList& GetDomainList() { return m_dom; }
	const FEDomainList& GetDomainList() const { return m_dom; }

	// Get an element
	FEElement& Element(int i);
	const FEElement& Element(int i) const;

	// create node list from this element set
	FENodeList GetNodeList() const;

public:
	void Serialize(DumpStream& ar);

	static void SaveClass(DumpStream& ar, FEElementSet* p);
	static FEElementSet* LoadClass(DumpStream& ar, FEElementSet* p);

private:
	// Build the lookup table
	void BuildLUT();

protected:
	std::vector<int>	m_Elem;		//!< list of elements' global ID

	FEDomainList		m_dom;	//!< domain list that generated the element set

	// used for fast lookup in GetLocalIndex
	std::vector<int>	m_LUT;
	int					m_minID, m_maxID;
};

inline int FEElementSet::GetLocalIndex(const FEElement& el) const
{
	int eid = el.GetID();
	if ((eid < m_minID) || (eid > m_maxID)) return -1;
	else return m_LUT[el.GetID() - m_minID];
}

inline bool FEElementSet::Contains(const FEElement& el) const
{
	return (GetLocalIndex(el) != -1);
}
