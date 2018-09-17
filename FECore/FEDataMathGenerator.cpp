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
	MObjBuilder mob;
	MSimpleExpression* val = dynamic_cast<MSimpleExpression*>(mob.Create(m_math, true));
	assert(val);
	if (val == 0) return false;

	MVariable* var_x = val->FindVariable("X");
	MVariable* var_y = val->FindVariable("Y");
	MVariable* var_z = val->FindVariable("Z");

	int N = set.size();
	ar.Create(N);
	for (int i = 0; i<N; ++i)
	{
		const FENode* ni = set.Node(i);

		vec3d ri = ni->m_r0;
		if (var_x) var_x->value(ri.x);
		if (var_y) var_y->value(ri.x);
		if (var_z) var_z->value(ri.x);

		double vi = val->value();
		ar.setValue(i, vi);
	}

	delete val;

	return true;
}

// generate the data array for the given facet set
bool FEDataMathGenerator::Generate(FESurfaceMap& data, const FEFacetSet& surf)
{
	MObjBuilder mob;
	MSimpleExpression* val = dynamic_cast<MSimpleExpression*>(mob.Create(m_math, true));
	assert(val);
	if (val == 0) return false;

	MVariable* var_x = val->FindVariable("X");
	MVariable* var_y = val->FindVariable("Y");
	MVariable* var_z = val->FindVariable("Z");

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
			if (var_x) var_x->value(ri.x);
			if (var_y) var_y->value(ri.x);
			if (var_z) var_z->value(ri.x);

			double vi = val->value();

			data.setValue(i, j, vi);
		}
	}

	return true;
}
