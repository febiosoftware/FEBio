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
#include <vector>
#include "fecore_api.h"

class FEMesh;
class FEElementList;
class FEElemElemList;
class FEEdgeList;

class FECORE_API FEFaceList
{
public:
	struct FACE
	{
		int		ntype;		// 3 = triangle, 4 = quad
		int		node[4];
		int		nsurf;		// 1 if facet is on surface 
		int		nbr[4];		// neighbor list

		bool IsEqual(int* n) const;

		bool HasEdge(int a, int b) const;
	};

public:
	FEFaceList();
	FEFaceList(const FEFaceList& faceList);

	bool Create(FEMesh& mesh, FEElemElemList& ENL);

	int Faces() const;

	const FACE& operator [] (int i) const;
	const FACE& Face(int i) const;

	FEMesh* GetMesh();

	// Extract the surface only
	FEFaceList GetSurface() const;

	// build the neighbor list
	void BuildNeighbors();

protected:
	FEMesh*				m_mesh;
	std::vector<FACE>	m_faceList;
};

class FECORE_API FENodeFaceList
{
public:
	FENodeFaceList();

	bool Create(FEFaceList& FL);

	int Faces(int node) const;
	const std::vector<int>& FaceList(int node) const;

private:
	std::vector<std::vector<int> >	m_NFL;
};

class FECORE_API FEElementFaceList
{
public:
	FEElementFaceList();

	bool Create(FEElementList& elemList, FEFaceList& faceList);

	int Faces(int elem) const;
	const std::vector<int>& FaceList(int elem) const;

private:
	std::vector<std::vector<int> >	m_EFL;
};

class FECORE_API FEFaceEdgeList
{
public:
	FEFaceEdgeList();

	bool Create(FEFaceList& faceList, FEEdgeList& edgeList);

	int Edges(int nface);
	const std::vector<int>& EdgeList(int nface) const;

private:
	std::vector<std::vector<int> > m_FEL;
};

class FECORE_API FENodeEdgeList
{
public:
	FENodeEdgeList();

	bool Create(FEEdgeList& edgeList);

	const std::vector<int>& EdgeList(int node) const { return m_NEL[node]; }

private:
	std::vector<std::vector<int> > m_NEL;
};
