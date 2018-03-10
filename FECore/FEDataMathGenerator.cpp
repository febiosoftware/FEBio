#include "stdafx.h"
#include "FEDataMathGenerator.h"
#include "MathParser.h"
#include "FEMesh.h"

//=======================================================================================
FEDataMathGenerator::FEDataMathGenerator()
{
}

// set the math expression
void FEDataMathGenerator::setExpression(const std::string& math)
{
	m_math = math;
}

// generate the data array for the given node set
bool FEDataMathGenerator::Generate(FENodeDataMap& ar, const FENodeSet& set)
{
	MathParser parser;

	int N = set.size();
	ar.Create(N);
	int ierr;
	for (int i = 0; i<N; ++i)
	{
		const FENode* ni = set.Node(i);

		vec3d ri = ni->m_r0;
		parser.SetVariable("X", ri.x);
		parser.SetVariable("Y", ri.y);
		parser.SetVariable("Z", ri.z);

		double vi = parser.eval(m_math.c_str(), ierr);
		if (ierr != 0) return false;

		ar.setValue(i, vi);
	}

	return true;
}

// generate the data array for the given facet set
bool FEDataMathGenerator::Generate(FESurfaceMap& data, const FEFacetSet& surf)
{
	MathParser parser;

	const FEMesh& mesh = *surf.GetMesh();

	int N = surf.Faces();
	data.Create(&surf);
	int ierr;
	for (int i = 0; i<N; ++i)
	{
		const FEFacetSet::FACET& face = surf.Face(i);

		int nf = face.ntype;
		for (int j=0; j<nf; ++j)
		{
			vec3d ri = mesh.Node(face.node[j]).m_r0;
			parser.SetVariable("X", ri.x);
			parser.SetVariable("Y", ri.y);
			parser.SetVariable("Z", ri.z);

			double vi = parser.eval(m_math.c_str(), ierr);
			if (ierr != 0) return false;

			data.setValue(i, j, vi);
		}
	}

	return true;
}
