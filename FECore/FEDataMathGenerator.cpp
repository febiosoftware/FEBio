#include "stdafx.h"
#include "FEDataMathGenerator.h"
#include "MathObject.h"
#include "MObjBuilder.h"
#include "FEMesh.h"

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
	MVariable* var_x = m_val.AddVariable("X");
	MVariable* var_y = m_val.AddVariable("Y");
	MVariable* var_z = m_val.AddVariable("Z");
	if (m_val.Create(m_math) == false) return false;

	return true;
}

double FEDataMathGenerator::value(const vec3d& r)
{
	vector<double> p{r.x, r.y, r.z};
	return m_val.value_s(p);
}
