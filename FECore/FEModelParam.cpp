#include "stdafx.h"
#include "FEModelParam.h"
#include "MObjBuilder.h"
#include "FEDataArray.h"

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

FEMappedValue::FEMappedValue(FEDomain* dom, FEDataMap* val)
{
	m_dom = dom;
	m_val = val;
}

double FEMappedValue::eval(const FEMaterialPoint& pt)
{
	return m_val->value(pt);
}

FEModelParam::FEModelParam() 
{ 
	m_scl = 1.0;
	m_val = new FEConstValue(0.0); 
}

// set the value
void FEModelParam::setValue(double v)
{ 
	setValuator(new FEConstValue(v));
}

// set the valuator
void FEModelParam::setValuator(FEValuator* val)
{
	if (m_val) delete m_val;
	m_val = val;
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
