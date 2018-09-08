#include "stdafx.h"
#include "FEMatParam.h"

FEMatExpression::FEMatExpression(const std::string& s) : m_expr(s)
{
}

double FEMatExpression::eval(const FEMaterialPoint& pt)
{
	m_math.SetVariable("X", pt.m_r0.x);
	m_math.SetVariable("Y", pt.m_r0.y);
	m_math.SetVariable("Z", pt.m_r0.z);

	const char* sz = m_expr.c_str();
	int ierr;
	return m_math.eval(sz, ierr);
}

FEMatMappedValue::FEMatMappedValue(FEDomain* dom, std::vector<double>& values)
{
	m_dom = dom;
	m_val = values;
}

double FEMatMappedValue::eval(const FEMaterialPoint& pt)
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

FEMaterialParam::FEMaterialParam() { m_val = new FEMatConstValue(0.0); }

// set the value
void FEMaterialParam::setValue(double v)
{ 
	setValuator(new FEMatConstValue(v));
}

// set the valuator
void FEMaterialParam::setValuator(FEMatValuator* val)
{
	if (m_val) delete m_val;
	m_val = val;
}
