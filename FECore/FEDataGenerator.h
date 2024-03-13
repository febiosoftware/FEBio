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
#include "vec3d.h"
#include "FEModelComponent.h"

class FENodeSet;
class FEEdgeList;
class FEFacetSet;
class FEElementSet;
class FEDataMap;

// Data generators are used to generate mesh data sections algorithmically
class FECORE_API FEMeshDataGenerator : public FEModelComponent
{
	FECORE_SUPER_CLASS(FEMESHDATAGENERATOR_ID)
	FECORE_BASE_CLASS(FEMeshDataGenerator)

public:
	FEMeshDataGenerator(FEModel* fem);
	virtual ~FEMeshDataGenerator();

	// this function gives the data generator a chance to initialize itself
	// and check for any input problems.
	virtual bool Init();

	// evaluate the data at specific time
	virtual void Evaluate(double time);

	// generate the mesh data section
	virtual FEDataMap* Generate() = 0;
};

// class for generating data on node sets
class FECORE_API FENodeDataGenerator : public FEMeshDataGenerator
{
	FECORE_BASE_CLASS(FENodeDataGenerator)

public:
	FENodeDataGenerator(FEModel* fem);

	void SetNodeSet(FENodeSet* nodeSet);

	FENodeSet* GetNodeSet();

protected:
	FENodeSet* m_nodeSet;
};

// class for generating data on edges
class FECORE_API FEEdgeDataGenerator : public FEMeshDataGenerator
{
	FECORE_BASE_CLASS(FEEdgeDataGenerator)

public:
	FEEdgeDataGenerator(FEModel* fem);

	void SetEdgeList(FEEdgeList* edgeSet);

	FEEdgeList* GetEdgeList();

protected:
	FEEdgeList* m_edgeList;
};

// class for generating data on surfaces
class FECORE_API FEFaceDataGenerator : public FEMeshDataGenerator
{
	FECORE_BASE_CLASS(FEFaceDataGenerator)

public:
	FEFaceDataGenerator(FEModel* fem);

	void SetFacetSet(FEFacetSet* surf);
	FEFacetSet* GetFacetSet();

private:
	FEFacetSet* m_surf;
};

// class for generating data on element sets
class FECORE_API FEElemDataGenerator : public FEMeshDataGenerator
{
	FECORE_BASE_CLASS(FEElemDataGenerator)

public:
	FEElemDataGenerator(FEModel* fem);

	void SetElementSet(FEElementSet* elset);
	FEElementSet* GetElementSet();

private:
	FEElementSet* m_elemSet;
};
