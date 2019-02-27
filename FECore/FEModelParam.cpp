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

FEModelParam::~FEModelParam()
{
}

// serialization
void FEModelParam::Serialize(DumpStream& ar)
{
	ar & m_scl;
}

//---------------------------------------------------------------------------------------
FEParamDouble::FEParamDouble()
{
	m_val = fecore_new<FEScalarValuator>("const", nullptr);
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
	FEConstValue* val = fecore_new<FEConstValue>("const", nullptr);
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

// is this a const value
bool FEParamDouble::isConst() const { return m_val->isConst(); };

// get the const value (returns 0 if param is not const)
double& FEParamDouble::constValue() { return *m_val->constValue(); }
double FEParamDouble::constValue() const { return *m_val->constValue(); }

void FEParamDouble::Serialize(DumpStream& ar)
{
	FEModelParam::Serialize(ar);
	ar & m_val;
}

//---------------------------------------------------------------------------------------
FEParamVec3::FEParamVec3()
{
	m_val = fecore_new<FEVec3dValuator>("vector", nullptr);
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
	FEConstValueVec3* val = fecore_new<FEConstValueVec3>("vector", nullptr);
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

void FEParamVec3::Serialize(DumpStream& ar)
{
	FEModelParam::Serialize(ar);
	ar & m_val;
}

//==========================================================================

FEParamMat3d::FEParamMat3d()
{
	m_val = fecore_new<FEMat3dValuator>("const", nullptr);
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
	FEConstValueMat3d* val = fecore_new<FEConstValueMat3d>("const", nullptr);
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

void FEParamMat3d::Serialize(DumpStream& ar)
{
	FEModelParam::Serialize(ar);
	ar & m_val;
}
