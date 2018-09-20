#include "stdafx.h"
#include "FEModelParam.h"
#include "MObjBuilder.h"
#include "FEDataArray.h"

//---------------------------------------------------------------------------------------
FEMathExpression::FEMathExpression(const std::string& s) : m_expr(s)
{
	m_math.AddVariable("X");
	m_math.AddVariable("Y");
	m_math.AddVariable("Z");
	bool b = m_math.Create(s);
	assert(b);
}

FEMathExpression::~FEMathExpression()
{
}

double FEMathExpression::eval(const FEMaterialPoint& pt)
{
	std::vector<double> var(3);
	var[0] = pt.m_r0.x;
	var[1] = pt.m_r0.y;
	var[2] = pt.m_r0.z;
	return m_math.value_s(var);
}

//---------------------------------------------------------------------------------------
FEMappedValue::FEMappedValue(FEDomain* dom, FEDataMap* val)
{
	m_dom = dom;
	m_val = val;
}

double FEMappedValue::eval(const FEMaterialPoint& pt)
{
	return m_val->value(pt);
}

//---------------------------------------------------------------------------------------
FEModelParam::FEModelParam()
{ 
	m_dom = 0;
	m_scl = 1.0;
}

// set the domain
void FEModelParam::setDomain(FEDomain* dom)
{
	m_dom = dom;
}

// get the domain
FEDomain* FEModelParam::getDomain()
{
	return m_dom;
}

//---------------------------------------------------------------------------------------
FEParamDouble::FEParamDouble()
{
	m_val = new FEConstValue(0.0);
}

// set the value
void FEParamDouble::setValue(double v)
{
	setValuator(new FEConstValue(v));
}

// set the valuator
void FEParamDouble::setValuator(FEValuator<double>* val)
{
	if (m_val) delete m_val;
	m_val = val;
}

//=======================================================================================

//---------------------------------------------------------------------------------------
FEMathExpressionVec3::FEMathExpressionVec3(const std::string& sx, const std::string& sy, const std::string& sz)
{
	for (int i = 0; i < 3; ++i)
	{
		m_math[i].AddVariable("X");
		m_math[i].AddVariable("Y");
		m_math[i].AddVariable("Z");
	}
	bool b;
	b = m_math[0].Create(sx); assert(b);
	b = m_math[1].Create(sy); assert(b);
	b = m_math[2].Create(sz); assert(b);
}

vec3d FEMathExpressionVec3::eval(const FEMaterialPoint& pt)
{
	std::vector<double> var(3);
	var[0] = pt.m_r0.x;
	var[1] = pt.m_r0.y;
	var[2] = pt.m_r0.z;
	double vx = m_math[0].value_s(var);
	double vy = m_math[1].value_s(var);
	double vz = m_math[2].value_s(var);
	return vec3d(vx, vy, vz);
}


//---------------------------------------------------------------------------------------
FEMappedValueVec3::FEMappedValueVec3(FEDomain* dom, FEDataMap* val)
{
	m_dom = dom;
	m_val = val;
}

vec3d FEMappedValueVec3::eval(const FEMaterialPoint& pt)
{
	return m_val->valueVec3d(pt);
}

//---------------------------------------------------------------------------------------
FEParamVec3::FEParamVec3()
{
	m_val = new FEConstValueVec3(vec3d(0,0,0));
}

// set the value
void FEParamVec3::setValue(const vec3d& v)
{
	setValuator(new FEConstValueVec3(v));
}

// set the valuator
void FEParamVec3::setValuator(FEValuator<vec3d>* val)
{
	if (m_val) delete m_val;
	m_val = val;
}
