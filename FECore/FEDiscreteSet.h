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
#include <vector>
#include <string>

//-----------------------------------------------------------------------------
// Forward declarations
class FEMesh;
class DumpStream;

//-----------------------------------------------------------------------------
//! Defines a discrete element set (i.e. node-pairs)
class FECORE_API FEDiscreteSet
{
public:
	struct NodePair
	{
		int	n0, n1;

		void Serialize(DumpStream& ar);
	};

public:
	FEDiscreteSet(FEMesh* pm);
	void create(int n);
	int size() const { return (int)m_pair.size(); }

	void add(int n0, int n1);

	void SetName(const std::string& name);
	const std::string& GetName() const;

	const NodePair& Element(int i) const { return m_pair[i]; }

	void Serialize(DumpStream& ar);

private:
	FEMesh*					m_pmesh;
	std::vector<NodePair>	m_pair;		//!< list of discrete elements
	std::string				m_name;
};
