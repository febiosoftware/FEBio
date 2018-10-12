#include "stdafx.h"
#include "FEDataMathGenerator.h"
#include "MathObject.h"
#include "MObjBuilder.h"
#include "FEMesh.h"
using namespace std;

BEGIN_FECORE_CLASS(FEDataMathGenerator, FEDataGenerator)
	ADD_PARAMETER(m_math, "math");
END_FECORE_CLASS();

FEDataMathGenerator::FEDataMathGenerator(FEModel* fem) : FEDataGenerator(fem)
{
}

// set the math expression
void FEDataMathGenerator::setExpression(const std::string& math)
{
	m_math = math;
}

bool FEDataMathGenerator::Init()
{
	string tmp = m_math;

	// split the math string at ','
	vector<string> strings;
	size_t pos = 0;
	while ((pos = tmp.find(',')) != string::npos)
	{
		string t = tmp.substr(0, pos);
		strings.push_back(t);
		tmp.erase(0, pos + 1);
	}
	strings.push_back(tmp);

	if ((strings.size() == 0) || (strings.size() > 3)) return false;

	for (size_t i = 0; i < strings.size(); ++i)
	{
		MVariable* var_x = m_val[i].AddVariable("X");
		MVariable* var_y = m_val[i].AddVariable("Y");
		MVariable* var_z = m_val[i].AddVariable("Z");
		MSimpleExpression val; 
		if (val.Create(strings[i]) == false) return false;
		m_val.push_back(val);
	}

	return true;
}

void FEDataMathGenerator::value(const vec3d& r, double& data)
{
	vector<double> p{r.x, r.y, r.z};
	assert(m_val.size() == 1);
	data = m_val[0].value_s(p);
}

void FEDataMathGenerator::value(const vec3d& r, vec3d& data)
{
	vector<double> p{ r.x, r.y, r.z };
	assert(m_val.size() <= 3);
	data.x = m_val[0].value_s(p);
	data.y = m_val[1].value_s(p);
	data.z = m_val[2].value_s(p);
}
