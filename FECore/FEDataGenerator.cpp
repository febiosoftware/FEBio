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
	for (int i = 0; i<N; ++i)
	{
		const FENode* ni = set.Node(i);
		vec3d ri = ni->m_r0;
		double vi = value(ri);
		ar.setValue(i, vi);
	}

	return true;
}

// generate the data array for the given facet set
bool FEDataGenerator::Generate(FESurfaceMap& data, const FEFacetSet& surf)
{
	int ntype = data.DataSize();
	const FEMesh& mesh = *surf.GetMesh();

	int N = surf.Faces();
	data.Create(&surf);
	for (int i = 0; i<N; ++i)
	{
		const FEFacetSet::FACET& face = surf.Face(i);

		int nf = face.ntype;
		for (int j=0; j<nf; ++j)
		{
			vec3d ri = mesh.Node(face.node[j]).m_r0;
			if (ntype == 1)
			{
				double vx = value(ri);
				data.setValue(i, j, vx);
			}
			else if (ntype == 2)
			{
			}
			else if (ntype == 3)
			{
			}
		}
	}
	return true;
}

// generate the data array for the given element set
bool FEDataGenerator::Generate(FEDomainMap& data, FEElementSet& set)
{
	int ntype = data.DataSize();
	FEMesh& mesh = *set.GetMesh();

	int N = set.Elements();
	data.Create(&set);
	for (int i = 0; i<N; ++i)
	{
		FEElement& el = *mesh.FindElementFromID(set[i]);

		int ne = el.Nodes();
		for (int j = 0; j < ne; ++j)
		{
			vec3d ri = mesh.Node(el.m_node[j]).m_r0;

			if (ntype == 1)
			{
				double vx = value(ri);
				data.setValue(i, j, vx);
			}
			else if (ntype == 2)
			{
			}
			else if (ntype == 3)
			{
			}
		}
	}

	return true;
}

