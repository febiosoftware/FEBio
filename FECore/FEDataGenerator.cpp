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



#include "stdafx.h"
#include "FEDataGenerator.h"
#include "FEMesh.h"
#include "FENodeDataMap.h"
#include "FEDomainMap.h"
#include "FEElementSet.h"
#include "log.h"

FEMeshDataGenerator::FEMeshDataGenerator(FEModel* fem) : FEModelComponent(fem)
{
}

FEMeshDataGenerator::~FEMeshDataGenerator()
{
}

bool FEMeshDataGenerator::Init()
{
	return true;
}

void FEMeshDataGenerator::Evaluate(double time)
{

}

//-----------------------------------------------------------------------------
FENodeDataGenerator::FENodeDataGenerator(FEModel* fem) : FEMeshDataGenerator(fem)
{
	m_nodeSet = nullptr;
}

void FENodeDataGenerator::SetNodeSet(FENodeSet* nodeSet)
{
	m_nodeSet = nodeSet;
}

FENodeSet* FENodeDataGenerator::GetNodeSet()
{
	return m_nodeSet;
}

//-----------------------------------------------------------------------------
FEEdgeDataGenerator::FEEdgeDataGenerator(FEModel* fem) : FEMeshDataGenerator(fem)
{
	m_edgeList = nullptr;
}

void FEEdgeDataGenerator::SetEdgeList(FEEdgeList* edgeSet)
{
	m_edgeList = edgeSet;
}

FEEdgeList* FEEdgeDataGenerator::GetEdgeList()
{
	return m_edgeList;
}

//-----------------------------------------------------------------------------
FEFaceDataGenerator::FEFaceDataGenerator(FEModel* fem) : FEMeshDataGenerator(fem)
{
	m_surf = nullptr;
}

void FEFaceDataGenerator::SetFacetSet(FEFacetSet* surf)
{
	m_surf = surf;
}

FEFacetSet* FEFaceDataGenerator::GetFacetSet()
{
	return m_surf;
}

//-----------------------------------------------------------------------------
FEElemDataGenerator::FEElemDataGenerator(FEModel* fem) : FEMeshDataGenerator(fem)
{
	m_elemSet = nullptr;
}

void FEElemDataGenerator::SetElementSet(FEElementSet* elset) { m_elemSet = elset; }
FEElementSet* FEElemDataGenerator::GetElementSet() { return m_elemSet; }
