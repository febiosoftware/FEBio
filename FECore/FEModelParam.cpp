#include "stdafx.h"
#include "FEModelParam.h"
#include "MObjBuilder.h"
#include "FEDataArray.h"

//---------------------------------------------------------------------------------------
FEMathExpression::FEMathExpression(const std::string& s, FECoreBase* pc) : m_expr(s)
{
	m_math.AddVariable("X");
	m_math.AddVariable("Y");
	m_math.AddVariable("Z");
	bool b = m_math.Create(s, true);

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
}

FEMathExpression::~FEMathExpression()
{
}

double FEMathExpression::eval(const FEMaterialPoint& pt)
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
			case FE_PARAM_INT          : var[3 + i] = (double) pi->value<int>(); break;
			case FE_PARAM_DOUBLE       : var[3 + i] = pi->value<double>(); break;
			case FE_PARAM_DOUBLE_MAPPED: var[3 + i] = pi->value<FEParamDouble>()(pt); break;
			}
		}
	}
	return m_math.value_s(var);
}

//---------------------------------------------------------------------------------------
FEMappedValue::FEMappedValue(FEDataMap* val)
{
	m_val = val;
}

double FEMappedValue::eval(const FEMaterialPoint& pt)
{
	return m_val->value(pt);
}

//---------------------------------------------------------------------------------------
FENodeMappedValue::FENodeMappedValue(FENodeDataMap* val)
{
	m_val = val;
}

double FENodeMappedValue::eval(const FEMaterialPoint& pt)
{
	return m_val->getValue(pt.m_index);
}

//---------------------------------------------------------------------------------------
FEModelParam::FEModelParam()
{ 
	m_scl = 1.0;
	m_dom = 0;
}

//---------------------------------------------------------------------------------------
FEParamDouble::FEParamDouble()
{
	m_val = new FEConstValue(0.0);
}

// set the value
void FEParamDouble::operator = (double v)
{
	setValuator(new FEConstValue(v));
}

// set the valuator
void FEParamDouble::setValuator(FEValuator<double>* val)
{
	if (m_val) delete m_val;
	m_val = val;
}

// is this a const value
bool FEParamDouble::isConst() const
{
	return (dynamic_cast<FEConstValue*>(m_val) != 0);
}

// get the const value (returns 0 if param is not const)
double FEParamDouble::constValue() const
{
	FEConstValue* cv = dynamic_cast<FEConstValue*>(m_val);
	if (cv) return cv->value(); else return 0.0;
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
FEMappedValueVec3::FEMappedValueVec3(FEDataMap* val)
{
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
void FEParamVec3::operator = (const vec3d& v)
{
	setValuator(new FEConstValueVec3(v));
}

// set the valuator
void FEParamVec3::setValuator(FEValuator<vec3d>* val)
{
	if (m_val) delete m_val;
	m_val = val;
}
