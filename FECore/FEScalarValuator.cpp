#include "stdafx.h"
#include "FEScalarValuator.h"
#include "FEMaterialPoint.h"
#include "FEModelParam.h"

void FEMathValue::setMathString(const std::string& s)
{
	m_expr = s;
}

bool FEMathValue::create(FECoreBase* pc)
{
	m_math.AddVariable("X");
	m_math.AddVariable("Y");
	m_math.AddVariable("Z");
	bool b = m_math.Create(m_expr, true);

	// lookup all the other variables.
	if (m_math.Variables() > 3)
	{
		assert(pc);
		for (int i = 3; i < m_math.Variables(); ++i)
		{
			MVariable* vari = m_math.Variable(i);

			ParamString ps(vari->Name().c_str());
			FEParam* p = pc->FindParameter(ps);
			assert(p);
			assert((p->type() == FE_PARAM_DOUBLE_MAPPED) || (p->type() == FE_PARAM_DOUBLE) || (p->type() == FE_PARAM_INT));

			m_vars.push_back(p);
		}
	}

	assert(b);
	return b;
}

FEMathValue::~FEMathValue()
{
}

FEScalarValuator* FEMathValue::copy()
{
	FEMathValue* newExpr = new FEMathValue(GetFEModel());
	newExpr->m_expr = m_expr;
	newExpr->m_math = m_math;
	newExpr->m_vars = m_vars;
	return newExpr;
}

double FEMathValue::operator()(const FEMaterialPoint& pt)
{
	std::vector<double> var(3 + m_vars.size());
	var[0] = pt.m_r0.x;
	var[1] = pt.m_r0.y;
	var[2] = pt.m_r0.z;
	if (m_vars.empty() == false)
	{
		for (int i = 0; i < (int)m_vars.size(); ++i)
		{
			FEParam* pi = m_vars[i];
			switch (pi->type())
			{
			case FE_PARAM_INT: var[3 + i] = (double)pi->value<int>(); break;
			case FE_PARAM_DOUBLE: var[3 + i] = pi->value<double>(); break;
			case FE_PARAM_DOUBLE_MAPPED: var[3 + i] = pi->value<FEParamDouble>()(pt); break;
			}
		}
	}
	return m_math.value_s(var);
}

//---------------------------------------------------------------------------------------

FEMappedValue::FEMappedValue(FEModel* fem) : FEScalarValuator(fem), m_val(nullptr), m_scale(0.0)
{
}

void FEMappedValue::setDataMap(FEDataMap* val, double scl)
{
	m_val = val;
	m_scale = scl;
}

double FEMappedValue::operator()(const FEMaterialPoint& pt)
{
	return m_scale*m_val->value(pt);
}

FEScalarValuator* FEMappedValue::copy()
{
	FEMappedValue* map = new FEMappedValue(GetFEModel());
	map->setDataMap(m_val, m_scale);
	return map;
}

//---------------------------------------------------------------------------------------

FENodeMappedValue::FENodeMappedValue(FEModel* fem) : FEScalarValuator(fem), m_val(nullptr), m_scale(0.0)
{

}

void FENodeMappedValue::setDataMap(FENodeDataMap* val, double scale)
{
	m_val = val;
	m_scale = scale;
}

double FENodeMappedValue::operator()(const FEMaterialPoint& pt)
{
	return m_scale*m_val->getValue(pt.m_index);
}

FEScalarValuator* FENodeMappedValue::copy()
{
	FENodeMappedValue* map = new FENodeMappedValue(GetFEModel());
	map->setDataMap(m_val, m_scale);
	return map;
}
