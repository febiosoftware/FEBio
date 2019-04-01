/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, Columbia University, and others.

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
#include "FEElement.h"
#include <vector>
#include <string>

//-----------------------------------------------------------------------------
//! This class defines a set of segments. This can be used in the creation of edges.
class FECORE_API FESegmentSet
{
public:
	struct SEGMENT
	{
		int	node[FEElement::MAX_NODES];
		int	ntype;	//	2=line2

		void Serialize(DumpStream& ar);
	};

public:
	FESegmentSet(FEMesh* pm);

	void SetName(const std::string& name);
	const std::string& GetName() const;

	void Create(int n);

	int Segments() { return (int)m_Seg.size(); }
	SEGMENT& Segment(int i);

	void Serialize(DumpStream& ar);

private:
	vector<SEGMENT>	m_Seg;
	std::string		m_name;
	FEMesh*			m_mesh;
};
