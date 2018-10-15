#include "stdafx.h"
#include "FEDataGenerator.h"
#include "FEMesh.h"
#include "FENodeDataMap.h"
#include "FESurfaceMap.h"
#include "FEDomainMap.h"
#include "FEElementSet.h"
#include "log.h"

FEDataGenerator::FEDataGenerator(FEModel* fem) : FECoreBase(fem, FEDATAGENERATOR_ID)
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
bool FEDataGenerator::Generate(FENodeDataMap& map, const FENodeSet& set)
{
	int N = set.size();
	map.Create(N);
	vector<double> p(3, 0.0);

	FEDataType dataType = map.DataType();
	for (int i = 0; i<N; ++i)
	{
		const FENode* ni = set.Node(i);
		vec3d ri = ni->m_r0;
		switch (dataType)
		{
		case FE_DOUBLE: { double d; value(ri, d); map.setValue(i, d); } break;
		case FE_VEC2D : { vec2d  d; value(ri, d); map.setValue(i, d); } break;
		case FE_VEC3D : { vec3d  d; value(ri, d); map.setValue(i, d); } break;
		case FE_MAT3D : { mat3d  d; value(ri, d); map.setValue(i, d); } break;
		}
	}

	return true;
}

// generate the data array for the given facet set
bool FEDataGenerator::Generate(FESurfaceMap& map, const FEFacetSet& surf)
{
	const FEMesh& mesh = *surf.GetMesh();

	FEDataType dataType = map.DataType();

	int N = surf.Faces();
	map.Create(&surf);
	for (int i = 0; i<N; ++i)
	{
		const FEFacetSet::FACET& face = surf.Face(i);

		int nf = face.ntype;
		for (int j=0; j<nf; ++j)
		{
			vec3d ri = mesh.Node(face.node[j]).m_r0;
			switch (dataType)
			{
			case FE_DOUBLE: { double d; value(ri, d); map.setValue(i, j, d); } break;
			case FE_VEC2D : { vec2d  d; value(ri, d); map.setValue(i, j, d); } break;
			case FE_VEC3D : { vec3d  d; value(ri, d); map.setValue(i, j, d); } break;
			case FE_MAT3D : { mat3d  d; value(ri, d); map.setValue(i, j, d); } break;
			}
		}
	}
	return true;
}

// generate the data array for the given element set
bool FEDataGenerator::Generate(FEDomainMap& map, FEElementSet& set)
{
	FEMesh& mesh = *set.GetMesh();

	FEDataType dataType = map.DataType();
	int N = set.Elements();
	map.Create(&set);
	for (int i = 0; i<N; ++i)
	{
		FEElement& el = *mesh.FindElementFromID(set[i]);

		int ne = el.Nodes();
		for (int j = 0; j < ne; ++j)
		{
			vec3d ri = mesh.Node(el.m_node[j]).m_r0;
			switch (dataType)
			{
			case FE_DOUBLE: { double d; value(ri, d); map.setValue(i, j, d); } break;
			case FE_VEC2D : { vec2d  d; value(ri, d); map.setValue(i, j, d); } break;
			case FE_VEC3D : { vec3d  d; value(ri, d); map.setValue(i, j, d); } break;
			case FE_MAT3D : { mat3d  d; value(ri, d); map.setValue(i, j, d); } break;
			}
		}
	}
	return true;
}
