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
#include "FEStepComponent.h"
#include <FECore/FENodeList.h>
#include <FECore/FENodeSet.h>
#include "fecore_api.h"

class FECORE_API FEPeriodicLinearConstraint : public FEStepComponent
{
	struct NodeSetPair
	{
		FENodeList primary;
		FENodeList secondary;
	};

public:
	FEPeriodicLinearConstraint(FEModel* fem);
	~FEPeriodicLinearConstraint();

	void AddNodeSetPair(const FENodeList& ms, const FENodeList& ss, bool push_back = true);

	void SetReferenceNode(int node) { m_refNode = node; }

	void ExcludeNodes(const FENodeSet& ps) { m_exclude = ps; }

public:
	// generate the linear constraints
	bool GenerateConstraints(FEModel* fem);

private:
	std::vector<NodeSetPair>	m_set;	// list of node set pairs
	FENodeSet	m_exclude;		// nodes to exclude
	int			m_refNode;		// reference node
};
