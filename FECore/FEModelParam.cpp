#include "stdafx.h"
#include "FEModelParam.h"
#include "MObjBuilder.h"

FEMathExpression::FEMathExpression(const std::string& s) : m_expr(s)
{
	MObjBuilder mob;
	m_math = dynamic_cast<MSimpleExpression*>(mob.Create(s, true));
	assert(m_math);

	m_var[0] = m_math->FindVariable("X");
	m_var[1] = m_math->FindVariable("Y");
	m_var[2] = m_math->FindVariable("Z");
}

FEMathExpression::~FEMathExpression()
{
	delete m_math;
}

double FEMathExpression::eval(const FEMaterialPoint& pt)
{
	// TODO: This is not thread safe! Need a better mechanism. For instance, pass variable values as parameter to value() function.
	if (m_var[0]) m_var[0]->value(pt.m_r0.x);
	if (m_var[1]) m_var[1]->value(pt.m_r0.y);
	if (m_var[2]) m_var[2]->value(pt.m_r0.z);

	return m_math->value();
}

FEMappedValue::FEMappedValue(FEDomain* dom, std::vector<double>& values)
{
	m_dom = dom;
	m_val = values;
}

double FEMappedValue::eval(const FEMaterialPoint& pt)
{
	// get the element this material point is in
	FEElement* pe = pt.m_elem;
	assert(pe);

	// make sure this element belongs to this domain
	assert(pe->GetDomain() == m_dom);

	// get its local ID
	int lid = pe->GetLocalID();

	assert((lid >= 0) && (lid < (int)m_val.size()));

	return m_val[lid];
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
