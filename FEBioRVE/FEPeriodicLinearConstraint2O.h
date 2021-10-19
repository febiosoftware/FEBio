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
#include <FECore/FEMesh.h>
#include "febiorve_api.h"
#include <vector>

class FEModel;

class FEBIORVE_API FEPeriodicLinearConstraint2O
{
	class NodeSetSet
	{
	public:
		NodeSetSet();
		NodeSetSet(const NodeSetSet& nss);
		void operator = (const NodeSetSet& nns);

	public:
		FENodeList	primary;
		FENodeList	secondary;
	};

public:
	FEPeriodicLinearConstraint2O();
	~FEPeriodicLinearConstraint2O();

	void AddNodeSetPair(const FENodeList& ms, const FENodeList& ss, bool push_back = true);

	bool GenerateConstraints(FEModel* fem);

private:
	int closestNode(FEMesh& mesh, const FENodeList& set, const vec3d& r);
	void addLinearConstraint(FEModel& fem, int parent, int child);

private:
	std::vector<NodeSetSet>	m_set;	// list of node set pairs
};
