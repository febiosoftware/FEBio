#include "stdafx.h"
#include "FEDataMathGenerator.h"
#include "MathObject.h"
#include "MObjBuilder.h"
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
	MSimpleExpression val;
	MVariable* var_x = val.AddVariable("X");
	MVariable* var_y = val.AddVariable("Y");
	MVariable* var_z = val.AddVariable("Z");
	if (val.Create(m_math) == false) return false;

	int N = set.size();
	ar.Create(N);
	for (int i = 0; i<N; ++i)
	{
		const FENode* ni = set.Node(i);

		vec3d ri = ni->m_r0;
		var_x->value(ri.x);
		var_y->value(ri.x);
		var_z->value(ri.x);

		double vi = val.value();
		ar.setValue(i, vi);
	}

	return true;
}

// generate the data array for the given facet set
bool FEDataMathGenerator::Generate(FESurfaceMap& data, const FEFacetSet& surf)
{
	MSimpleExpression val;
	MVariable* var_x = val.AddVariable("X");
	MVariable* var_y = val.AddVariable("Y");
	MVariable* var_z = val.AddVariable("Z");
	if (val.Create(m_math) == false) return false;

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
			var_x->value(ri.x);
			var_y->value(ri.y);
			var_z->value(ri.z);

			double vi = val.value();

			data.setValue(i, j, vi);
		}
	}

	return true;
}
