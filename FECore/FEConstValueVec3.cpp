/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2020 University of Utah, The Trustees of Columbia University in
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
#include "FEConstValueVec3.h"
#include "FEMaterialPoint.h"
#include "FEMeshPartition.h"
#include "FENode.h"
#include "quatd.h"

//==================================================================================
BEGIN_FECORE_CLASS(FEConstValueVec3, FEVec3dValuator)
	ADD_PARAMETER(m_val, "vector");
END_FECORE_CLASS();

FEConstValueVec3::FEConstValueVec3(FEModel* fem) : FEVec3dValuator(fem) {}

FEVec3dValuator* FEConstValueVec3::copy()
{
	FEConstValueVec3* val = new FEConstValueVec3(GetFEModel());
	val->m_val = m_val;
	return val;
}

//==================================================================================
BEGIN_FECORE_CLASS(FEMathValueVec3, FEVec3dValuator)
	ADD_PARAMETER(m_expr, "math");
END_FECORE_CLASS();

FEMathValueVec3::FEMathValueVec3(FEModel* fem) : FEVec3dValuator(fem)
{
}

//---------------------------------------------------------------------------------------
bool FEMathValueVec3::Init()
{
	size_t c1 = m_expr.find(',', 0); if (c1 == string::npos) return false;
	size_t c2 = m_expr.find(',', c1 + 1); if (c2 == string::npos) return false;

	string sx = m_expr.substr(0, c1);
	string sy = m_expr.substr(c1 + 1, c2 - c1);
	string sz = m_expr.substr(c2 + 1, string::npos);

	return create(sx, sy, sz);
}

//---------------------------------------------------------------------------------------
bool FEMathValueVec3::create(const std::string& sx, const std::string& sy, const std::string& sz)
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

	return true;
}

vec3d FEMathValueVec3::operator()(const FEMaterialPoint& pt)
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
FEVec3dValuator* FEMathValueVec3::copy()
{
	FEMathValueVec3* newVal = new FEMathValueVec3(GetFEModel());
	newVal->m_math[0] = m_math[0];
	newVal->m_math[1] = m_math[1];
	newVal->m_math[2] = m_math[2];
	return newVal;
}

//---------------------------------------------------------------------------------------

FEMappedValueVec3::FEMappedValueVec3(FEModel* fem) : FEVec3dValuator(fem)
{
	m_val = nullptr;
}

void FEMappedValueVec3::setDataMap(FEDataMap* val, vec3d scl)
{
	m_val = val;
}

vec3d FEMappedValueVec3::operator()(const FEMaterialPoint& pt)
{
	vec3d r = m_val->valueVec3d(pt);
	return vec3d(r.x, r.y, r.z);
}

FEVec3dValuator* FEMappedValueVec3::copy()
{
	FEMappedValueVec3* map = new FEMappedValueVec3(GetFEModel());
	map->m_val = m_val;
	return map;
}

void FEMappedValueVec3::Serialize(DumpStream& ar)
{
	FEVec3dValuator::Serialize(ar);
	if (ar.IsShallow()) return;
	ar & m_val;
}

//=================================================================================================
BEGIN_FECORE_CLASS(FELocalVectorGenerator, FEVec3dValuator)
	ADD_PARAMETER(m_n, 2, "local");
END_FECORE_CLASS();

FELocalVectorGenerator::FELocalVectorGenerator(FEModel* fem) : FEVec3dValuator(fem)
{
	m_n[0] = m_n[1] = 0;
}

bool FELocalVectorGenerator::Init()
{
	if ((m_n[0] == 0) && (m_n[1] == 0)) m_n[1] = 1;
	return FEVec3dValuator::Init();
}

vec3d FELocalVectorGenerator::operator () (const FEMaterialPoint& mp)
{
	FEElement* el = mp.m_elem; assert(el);

	FEMeshPartition* dom = el->GetMeshPartition();
	vec3d r0 = dom->Node(el->m_lnode[m_n[0]]).m_r0;
	vec3d r1 = dom->Node(el->m_lnode[m_n[1]]).m_r0;

	vec3d n = r1 - r0;
	n.unit();

	return n;
}

FEVec3dValuator* FELocalVectorGenerator::copy()
{
	FELocalVectorGenerator* map = new FELocalVectorGenerator(GetFEModel());
	map->m_n[0] = m_n[0];
	map->m_n[1] = m_n[1];
	return map;
}

//=================================================================================================
BEGIN_FECORE_CLASS(FESphericalVectorGenerator, FEVec3dValuator)
	ADD_PARAMETER(m_center, "center");
	ADD_PARAMETER(m_vector, "vector");
END_FECORE_CLASS();

FESphericalVectorGenerator::FESphericalVectorGenerator(FEModel* fem) : FEVec3dValuator(fem)
{
	m_center = vec3d(0, 0, 0);
	m_vector = vec3d(1, 0, 0);
}

bool FESphericalVectorGenerator::Init()
{
	// Make sure the vector is a unit vector
	m_vector.unit();
	return true;
}

FEVec3dValuator* FESphericalVectorGenerator::copy()
{
	FESphericalVectorGenerator* map = new FESphericalVectorGenerator(GetFEModel());
	map->m_center = m_center;
	map->m_vector = m_vector;
	return map;
}

vec3d FESphericalVectorGenerator::operator () (const FEMaterialPoint& mp)
{
	vec3d a = mp.m_r0 - m_center;
	a.unit();

	// setup the rotation
	vec3d e1(1, 0, 0);
	quatd q(e1, a);

	vec3d v = m_vector;
	//	v.unit();	
	q.RotateVector(v);

	return v;
}

//=================================================================================================
BEGIN_FECORE_CLASS(FECylindricalVectorGenerator, FEVec3dValuator)
	ADD_PARAMETER(m_center, "center");
	ADD_PARAMETER(m_axis, "axis")
	ADD_PARAMETER(m_vector, "vector");
END_FECORE_CLASS();

FECylindricalVectorGenerator::FECylindricalVectorGenerator(FEModel* fem) : FEVec3dValuator(fem)
{
	m_center = vec3d(0, 0, 0);
	m_axis = vec3d(0, 0, 1);
	m_vector = vec3d(1, 0, 0);
}

bool FECylindricalVectorGenerator::Init()
{
	// Make sure the axis and vector are unit vectors
	m_axis.unit();
	m_vector.unit();
	return true;
}

vec3d FECylindricalVectorGenerator::operator () (const FEMaterialPoint& mp)
{
	vec3d p = mp.m_r0 - m_center;

	// find the vector to the axis
	vec3d b = p - m_axis * (m_axis*p);
	b.unit();

	// setup the rotation
	vec3d e1(1, 0, 0);
	quatd q(e1, b);

	vec3d r = m_vector;
	//	r.unit();	
	q.RotateVector(r);

	return r;
}

FEVec3dValuator* FECylindricalVectorGenerator::copy()
{
	FECylindricalVectorGenerator* map = new FECylindricalVectorGenerator(GetFEModel());
	map->m_center = m_center;
	map->m_axis = m_axis;
	map->m_vector = m_vector;
	return map;
}


//=================================================================================================
BEGIN_FECORE_CLASS(FESphericalAnglesVectorGenerator, FEVec3dValuator)
	ADD_PARAMETER(m_theta, "theta");
	ADD_PARAMETER(m_phi, "phi");
END_FECORE_CLASS();

FESphericalAnglesVectorGenerator::FESphericalAnglesVectorGenerator(FEModel* fem) : FEVec3dValuator(fem)
{
	// equal to x-axis (1,0,0)
	m_theta = 0.0;
	m_phi = 90.0;
}

vec3d FESphericalAnglesVectorGenerator::operator () (const FEMaterialPoint& mp)
{
	// convert from degress to radians
	const double the = m_theta(mp)* PI / 180.;
	const double phi = m_phi(mp)* PI / 180.;

	// the fiber vector
	vec3d a;
	a.x = cos(the)*sin(phi);
	a.y = sin(the)*sin(phi);
	a.z = cos(phi);

	return a;
}

FEVec3dValuator* FESphericalAnglesVectorGenerator::copy()
{
	FESphericalAnglesVectorGenerator* v = fecore_alloc(FESphericalAnglesVectorGenerator, GetFEModel());
	v->m_theta = m_theta;
	v->m_phi = m_phi;
	return v;
}
