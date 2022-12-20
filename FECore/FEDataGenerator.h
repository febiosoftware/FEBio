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

//-----------------------------------------------------------------------------
class FENodeSet;
class FEFacetSet;
class FEElementSet;
class FENodeDataMap;
class FESurfaceMap;
class FEDomainMap;

//-----------------------------------------------------------------------------
// Data generators are used to generate values of model parameters. 
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

public:
	// overload  one of these functions for custom generators
	virtual void value(const vec3d& r, double& data) {}
	virtual void value(const vec3d& r, vec2d& data) {}
	virtual void value(const vec3d& r, vec3d& data) {}
	virtual void value(const vec3d& r, mat3d& data) {}
    virtual void value(const vec3d& r, mat3ds& data) {}
};

//-----------------------------------------------------------------------------
class FECORE_API FENodeDataGenerator : public FEMeshDataGenerator
{
	FECORE_BASE_CLASS(FENodeDataGenerator)

public:
	FENodeDataGenerator(FEModel* fem);

	// generate the data array for the given node set
	virtual bool Generate(FENodeDataMap& ar);

	virtual FENodeDataMap* Generate();

public:
	void SetNodeSet(FENodeSet* nodeSet);

	FENodeSet* GetNodeSet();

protected:
	FENodeSet* m_nodeSet;
};

//-----------------------------------------------------------------------------
class FECORE_API FEEdgeDataGenerator : public FEMeshDataGenerator
{
	FECORE_BASE_CLASS(FEEdgeDataGenerator)

public:
	FEEdgeDataGenerator(FEModel* fem);

	// generate the data array for the given node set
//	virtual bool Generate(FEEdgeDataMap& ar);
};

//-----------------------------------------------------------------------------
class FECORE_API FEFaceDataGenerator : public FEMeshDataGenerator
{
	FECORE_BASE_CLASS(FEFaceDataGenerator)

public:
	FEFaceDataGenerator(FEModel* fem);

	// generate the data array for the given facet set
	virtual bool Generate(FESurfaceMap& data);

	virtual FESurfaceMap* Generate();

public:
	void SetFacetSet(FEFacetSet* surf);
	FEFacetSet* GetFacetSet();

private:
	FEFacetSet* m_surf;
};

//-----------------------------------------------------------------------------
class FECORE_API FEElemDataGenerator : public FEMeshDataGenerator
{
	FECORE_BASE_CLASS(FEElemDataGenerator)

public:
	FEElemDataGenerator(FEModel* fem);

	// generate the data array for the given element set
	virtual bool Generate(FEDomainMap& data);

	virtual FEDomainMap* Generate();

public:
	void SetElementSet(FEElementSet* elset);
	FEElementSet* GetElementSet();

private:
	FEElementSet* m_elemSet;
};
