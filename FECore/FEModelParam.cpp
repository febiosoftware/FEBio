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
#include "FEModelParam.h"
#include "MObjBuilder.h"
#include "FEDataArray.h"
#include "DumpStream.h"
#include "FEConstValueVec3.h"

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
	if (ar.IsShallow() == false)
	{
		ar & m_dom;
	}
}

//---------------------------------------------------------------------------------------
FEParamDouble::FEParamDouble()
{
	m_val = fecore_new<FEScalarValuator>("const", nullptr);
}

FEParamDouble::~FEParamDouble()
{
	delete m_val;
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

void FEParamDouble::operator = (const FEParamDouble& p)
{
	if (m_val) delete m_val;
	m_val = p.m_val->copy();
	m_scl = p.m_scl;
//	m_dom = p.m_dom;
}

// set the valuator
void FEParamDouble::setValuator(FEScalarValuator* val)
{
	if (m_val) delete m_val;
	m_val = val;
	if (val) val->SetModelParam(this);
}

// get the valuator
FEScalarValuator* FEParamDouble::valuator()
{
	return m_val;
}

// is this a const value
bool FEParamDouble::isConst() const { return m_val->isConst(); };

// get the const value (returns 0 if param is not const)
double& FEParamDouble::constValue() { assert(isConst());  return *m_val->constValue(); }
double FEParamDouble::constValue() const { assert(isConst()); return *m_val->constValue(); }

void FEParamDouble::Serialize(DumpStream& ar)
{
	FEModelParam::Serialize(ar);
	ar & m_val;
}

bool FEParamDouble::Init()
{
	return (m_val ? m_val->Init() : true);
}

//---------------------------------------------------------------------------------------
FEParamVec3::FEParamVec3()
{
	m_val = fecore_new<FEVec3dValuator>("vector", nullptr);
}

FEParamVec3::~FEParamVec3()
{
	delete m_val;
}

FEParamVec3::FEParamVec3(const FEParamVec3& p)
{
	m_val = p.m_val->copy();
	m_scl = p.m_scl;
	m_dom = p.m_dom;
}

bool FEParamVec3::Init()
{
	return (m_val ? m_val->Init() : true);
}

// set the value
void FEParamVec3::operator = (const vec3d& v)
{
	FEConstValueVec3* val = fecore_new<FEConstValueVec3>("vector", nullptr);
	val->value() = v;
	setValuator(val);
}

void FEParamVec3::operator = (const FEParamVec3& p)
{
	m_val = p.m_val->copy();
	m_scl = p.m_scl;
//	m_dom = p.m_dom;
}

// set the valuator
void FEParamVec3::setValuator(FEVec3dValuator* val)
{
	if (m_val) delete m_val;
	m_val = val;
	if (val) val->SetModelParam(this);
}

FEVec3dValuator* FEParamVec3::valuator()
{
	return m_val;
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

FEParamMat3d::~FEParamMat3d()
{
	delete m_val;
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

void FEParamMat3d::operator = (const FEParamMat3d& p)
{
	m_val = p.m_val->copy();
	m_scl = p.m_scl;
//	m_dom = p.m_dom;
}

bool FEParamMat3d::Init()
{
	return (m_val ? m_val->Init() : true);
}

// set the valuator
void FEParamMat3d::setValuator(FEMat3dValuator* val)
{
	if (m_val) delete m_val;
	m_val = val;
	if (val) val->SetModelParam(this);
}

// get the valuator
FEMat3dValuator* FEParamMat3d::valuator()
{
	return m_val;
}

void FEParamMat3d::Serialize(DumpStream& ar)
{
	FEModelParam::Serialize(ar);
	ar & m_val;
}

//==========================================================================
FEParamMat3ds::FEParamMat3ds()
{
	m_val = fecore_new<FEMat3dsValuator>("const", nullptr);
}

FEParamMat3ds::~FEParamMat3ds()
{
	delete m_val;
}

FEParamMat3ds::FEParamMat3ds(const FEParamMat3ds& p)
{
	m_val = p.m_val->copy();
	m_scl = p.m_scl;
	m_dom = p.m_dom;
}

// set the value
void FEParamMat3ds::operator = (const mat3ds& v)
{
	FEConstValueMat3ds* val = fecore_new<FEConstValueMat3ds>("const", nullptr);
	val->value() = v;
	setValuator(val);
}

void FEParamMat3ds::operator = (const FEParamMat3ds& p)
{
	m_val = p.m_val->copy();
	m_scl = p.m_scl;
//	m_dom = p.m_dom;
}

// set the valuator
void FEParamMat3ds::setValuator(FEMat3dsValuator* val)
{
	if (m_val) delete m_val;
	m_val = val;
	if (val) val->SetModelParam(this);
}

FEMat3dsValuator* FEParamMat3ds::valuator()
{
	return m_val;
}

void FEParamMat3ds::Serialize(DumpStream& ar)
{
	FEModelParam::Serialize(ar);
	ar & m_val;
}
