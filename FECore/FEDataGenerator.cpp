#include "stdafx.h"
#include "FEDataGenerator.h"
#include "FEMesh.h"
#include "FENodeDataMap.h"
#include "FESurfaceMap.h"
#include "FEDomainMap.h"
#include "FEElementSet.h"
#include "log.h"

FEDataGenerator::FEDataGenerator(FEModel* fem) : FECoreBase(FEDATAGENERATOR_ID), m_fem(fem)
{
}

FEDataGenerator::~FEDataGenerator()
{
}

bool FEDataGenerator::Init()
{
	return true;
}

// generate the data array for the given node set
bool FEDataGenerator::Generate(FENodeDataMap& ar, const FENodeSet& set)
{
	int N = set.size();
	ar.Create(N);
	vector<double> p(3, 0.0);

	int dataSize = ar.DataSize();
	vector<double> data(dataSize);

	for (int i = 0; i<N; ++i)
	{
		const FENode* ni = set.Node(i);
		vec3d ri = ni->m_r0;
		value(ri, data);
		switch (dataSize)
		{
		case 1: ar.setValue(i, data[0]); break;
		case 2: ar.setValue(i, vec2d(data[0], data[1])); break;
		case 3: ar.setValue(i, vec3d(data[0], data[1], data[2])); break;
		}
	}

	return true;
}

// generate the data array for the given facet set
bool FEDataGenerator::Generate(FESurfaceMap& map, const FEFacetSet& surf)
{
	const FEMesh& mesh = *surf.GetMesh();

	int dataSize = map.DataSize();
	vector<double> data(dataSize);

	int N = surf.Faces();
	map.Create(&surf);
	for (int i = 0; i<N; ++i)
	{
		const FEFacetSet::FACET& face = surf.Face(i);

		int nf = face.ntype;
		for (int j=0; j<nf; ++j)
		{
			vec3d ri = mesh.Node(face.node[j]).m_r0;
			value(ri, data);
			switch (dataSize)
			{
			case 1: map.setValue(i, j, data[0]); break;
			case 2: map.setValue(i, j, vec2d(data[0], data[1])); break;
			case 3: map.setValue(i, j, vec3d(data[0], data[1], data[2])); break;
			}
		}
	}
	return true;
}

// generate the data array for the given element set
bool FEDataGenerator::Generate(FEDomainMap& map, FEElementSet& set)
{
	FEMesh& mesh = *set.GetMesh();

	int dataSize = map.DataSize();
	vector<double> data(dataSize);

	int N = set.Elements();
	map.Create(&set);
	for (int i = 0; i<N; ++i)
	{
		FEElement& el = *mesh.FindElementFromID(set[i]);

		int ne = el.Nodes();
		for (int j = 0; j < ne; ++j)
		{
			vec3d ri = mesh.Node(el.m_node[j]).m_r0;
			value(ri, data);
			switch (dataSize)
			{
			case 1: map.setValue(i, j, data[0]); break;
			case 2: map.setValue(i, j, vec2d(data[0], data[1])); break;
			case 3: map.setValue(i, j, vec3d(data[0], data[1], data[2])); break;
			}
		}
	}
	return true;
}
