/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2020 University of Utah, The Trustees of Columbia University in
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
#include "FECoreBase.h"

class FEElementSet;
class FEElement;

//-----------------------------------------------------------------------------
class FECORE_API FEMeshAdaptorSelection
{
public:
	struct Item {
		int		m_elementIndex;
		double	m_scaleFactor;
	};

public:
	FEMeshAdaptorSelection() {}
	FEMeshAdaptorSelection(size_t size) : m_itemList(size) {}

	void resize(size_t newSize) { m_itemList.resize(newSize); }
	Item& operator [] (size_t item) { return m_itemList[item]; }
	const Item& operator [] (size_t item) const { return m_itemList[item]; }
	bool empty() const { return m_itemList.empty(); }
	size_t size() const { return m_itemList.size(); }
	void push_back(int elemIndex, double scale) { m_itemList.push_back(Item{ elemIndex, scale }); }

private:
	std::vector<Item>	m_itemList;
};

//-----------------------------------------------------------------------------
// This class is a helper class for use in the mesh adaptors. Its purpose is to select
// elements based on some criterion. This element list is then usually passed to the 
// mesh adaptor.
class FECORE_API FEMeshAdaptorCriterion : public FECoreBase
{
	FECORE_SUPER_CLASS

public:
	FEMeshAdaptorCriterion(FEModel* fem);

	void SetSort(bool b);

	void SetMaxElements(int m);

public:

	// return a list of elements that satisfy the criterion
	// The elements will be selected from the element set. If nullptr is passed
	// for the element set, the entire mesh will be processed
	virtual FEMeshAdaptorSelection GetElementSelection(FEElementSet* elset);

	// This function needs to be overridden in order to select some elements
	// that satisfy the selection criterion
	// return true if the element satisfies the criterion, otherwise false
	// If this function returns true, the elemVal parameter should be set
	// This is used to sort the element list
	virtual bool Check(FEElement& el, double& elemVal);

private:
	bool	m_sortList;		// sort the list
	int		m_maxelem;		// the max nr of elements to return (or 0 if don't care)

	DECLARE_FECORE_CLASS();
};
