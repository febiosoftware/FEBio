#include "stdafx.h"
#include "FEModelParam.h"
#include "MObjBuilder.h"
#include "FEDataArray.h"

//---------------------------------------------------------------------------------------
FEModelParam::FEModelParam()
{ 
	m_scl = 1.0;
	m_dom = 0;
}

//---------------------------------------------------------------------------------------
FEParamDouble::FEParamDouble()
{
	m_val = new FEConstValue(nullptr);
}

FEParamDouble::FEParamDouble(const FEParamDouble& p)
{
	m_val = p.m_val->copy();
	m_scl = p.m_scl;
	m_dom = p.m_dom;
}

// set the value
void FEParamDouble::operator = (double v)
{
	FEConstValue* val = new FEConstValue(nullptr);
	*val->constValue() = v;
	setValuator(val);
}

// set the valuator
void FEParamDouble::setValuator(FEScalarValuator* val)
{
	if (m_val) delete m_val;
	m_val = val;
	if (val) val->SetModelParam(this);
}

//---------------------------------------------------------------------------------------
FEParamVec3::FEParamVec3()
{
	m_val = new FEConstValueVec3(nullptr);
}

FEParamVec3::FEParamVec3(const FEParamVec3& p)
{
	m_val = p.m_val->copy();
	m_scl = p.m_scl;
	m_dom = p.m_dom;
}

// set the value
void FEParamVec3::operator = (const vec3d& v)
{
	FEConstValueVec3* val = new FEConstValueVec3(nullptr);
	val->value() = v;
	setValuator(val);
}

// set the valuator
void FEParamVec3::setValuator(FEVec3dValuator* val)
{
	if (m_val) delete m_val;
	m_val = val;
	if (val) val->SetModelParam(this);
}

//==========================================================================

FEParamMat3d::FEParamMat3d()
{
	m_val = new FEConstValueMat3d(nullptr);
}

FEParamMat3d::FEParamMat3d(const FEParamMat3d& p)
{
	m_val = p.m_val->copy();
	m_scl = p.m_scl;
	m_dom = p.m_dom;
}

// set the value
void FEParamMat3d::operator = (const mat3d& v)
{
	FEConstValueMat3d* val = new FEConstValueMat3d(nullptr);
	val->value() = v;
	setValuator(val);
}

// set the valuator
void FEParamMat3d::setValuator(FEMat3dValuator* val)
{
	if (m_val) delete m_val;
	m_val = val;
	if (val) val->SetModelParam(this);
}
