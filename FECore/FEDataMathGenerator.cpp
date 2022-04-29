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
#include "FEDataMathGenerator.h"
#include "FENodeDataMap.h"
#include "MathObject.h"
#include "MObjBuilder.h"
#include "FEMesh.h"
using namespace std;

BEGIN_FECORE_CLASS(FEDataMathGenerator, FENodeDataGenerator)
	ADD_PARAMETER(m_math, "math");
END_FECORE_CLASS();

FEDataMathGenerator::FEDataMathGenerator(FEModel* fem) : FENodeDataGenerator(fem)
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
		MSimpleExpression val; 
		MVariable* var_x = val.AddVariable("X");
		MVariable* var_y = val.AddVariable("Y");
		MVariable* var_z = val.AddVariable("Z");
		if (val.Create(strings[i]) == false) return false;
		m_val.push_back(val);
	}

	return true;
}

FENodeDataMap* FEDataMathGenerator::Generate()
{
	assert(m_nodeSet);
	if (m_nodeSet == nullptr) return nullptr;
	FENodeDataMap* map = new FENodeDataMap(FE_DOUBLE);
	map->Create(m_nodeSet);
	if (FENodeDataGenerator::Generate(*map) == false)
	{
		delete map;
		map = nullptr;
		assert(false);
	}
	return map;
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
