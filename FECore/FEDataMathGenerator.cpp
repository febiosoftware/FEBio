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
		if (m_val[i].Create(strings[i]) == false) return false;
	}

	return true;
}

void FEDataMathGenerator::value(const vec3d& r, vector<double>& data)
{
	vector<double> p{r.x, r.y, r.z};
	assert(data.size() <= 3);
	for (size_t i = 0; i < data.size(); ++i)
	{
		data[i] = m_val[i].value_s(p);
	}
}
