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
	m_math.clear();
	m_math.push_back(math);
}

void FEDataMathGenerator::setExpression(const std::vector<std::string>& math)
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
	if (val.Create(m_math[0]) == false) return false;

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
	MSimpleExpression val[3];
	int ntype = data.DataSize();
	if (ntype != m_math.size()) return false;
	if (ntype > 3) return false;
	for (int i = 0; i < ntype; ++i)
	{
		MVariable* var_x = val[i].AddVariable("X");
		MVariable* var_y = val[i].AddVariable("Y");
		MVariable* var_z = val[i].AddVariable("Z");
		if (val[i].Create(m_math[i]) == false) return false;
	}

	const FEMesh& mesh = *surf.GetMesh();

	int N = surf.Faces();
	data.Create(&surf);
	vector<double> p(3, 0.0);
	for (int i = 0; i<N; ++i)
	{
		const FEFacetSet::FACET& face = surf.Face(i);

		int nf = face.ntype;
		for (int j=0; j<nf; ++j)
		{
			vec3d ri = mesh.Node(face.node[j]).m_r0;
			p[0] = ri.x;
			p[1] = ri.y;
			p[2] = ri.z;

			if (ntype == 1)
			{
				double vx = val[0].value_s(p);
				data.setValue(i, j, vx);
			}
			else if (ntype == 2)
			{
				double vx = val[0].value_s(p);
				double vy = val[1].value_s(p);
				data.setValue(i, j, vec2d(vx, vy));
			}
			else if (ntype == 3)
			{
				double vx = val[0].value_s(p);
				double vy = val[1].value_s(p);
				double vz = val[2].value_s(p);
				data.setValue(i, j, vec3d(vx, vy, vz));
			}
		}
	}

	return true;
}

// generate the data array for the given element set
bool FEDataMathGenerator::Generate(FEDomainMap& data, FEElementSet& set)
{
	MSimpleExpression val[3];
	int ntype = data.DataSize();
	if (ntype != m_math.size()) return false;
	if (ntype > 3) return false;
	for (int i = 0; i < ntype; ++i)
	{
		MVariable* var_x = val[i].AddVariable("X");
		MVariable* var_y = val[i].AddVariable("Y");
		MVariable* var_z = val[i].AddVariable("Z");
		if (val[i].Create(m_math[i]) == false) return false;
	}

	FEMesh& mesh = *set.GetMesh();

	int N = set.Elements();
	data.Create(&set);
	vector<double> p(3, 0.0);
	for (int i = 0; i<N; ++i)
	{
		FEElement& el = *mesh.FindElementFromID(set[i]);

		int ne = el.Nodes();
		for (int j = 0; j < ne; ++j)
		{
			vec3d ri = mesh.Node(el.m_node[j]).m_r0;
			p[0] = ri.x;
			p[1] = ri.y;
			p[2] = ri.z;

			if (ntype == 1)
			{
				double vx = val[0].value_s(p);
				data.setValue(i, j, vx);
			}
			else if (ntype == 2)
			{
				double vx = val[0].value_s(p);
				double vy = val[1].value_s(p);
				data.setValue(i, j, vec2d(vx, vy));
			}
			else if (ntype == 3)
			{
				double vx = val[0].value_s(p);
				double vy = val[1].value_s(p);
				double vz = val[2].value_s(p);
				data.setValue(i, j, vec3d(vx, vy, vz));
			}
		}
	}

	return true;
}
