#include "stdafx.h"
#include "FEMat3dValuator.h"
#include "FEModel.h"
#include "FEMesh.h"
#include "FEDataMap.h"

//=============================================================================
// FELocalMap
//-----------------------------------------------------------------------------

BEGIN_FECORE_CLASS(FEConstValueMat3d, FEMat3dValuator)
	ADD_PARAMETER(m_val, "const");
END_FECORE_CLASS();

FEConstValueMat3d::FEConstValueMat3d(FEModel* fem) : FEMat3dValuator(fem)
{
	m_val.zero();
}

FEMat3dValuator* FEConstValueMat3d::copy()
{
	FEConstValueMat3d* map = new FEConstValueMat3d(GetFEModel());
	map->m_val = m_val;
	return map;
}

//=============================================================================
// FELocalMap
//-----------------------------------------------------------------------------

BEGIN_FECORE_CLASS(FEMat3dLocalElementMap, FEMat3dValuator)
	ADD_PARAMETER(m_n, 3, "local");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEMat3dLocalElementMap::FEMat3dLocalElementMap(FEModel* pfem) : FEMat3dValuator(pfem)
{
	m_n[0] = 0;
	m_n[1] = 0;
	m_n[2] = 0;
}

//-----------------------------------------------------------------------------
bool FEMat3dLocalElementMap::Init()
{
	// check values
	if ((m_n[0] <= 0) && (m_n[1] <= 0) && (m_n[2] <= 0)) { m_n[0] = 1; m_n[1] = 2; m_n[2] = 4; }
	if (m_n[2] <= 0) m_n[2] = m_n[1];

	return true;
}

//-----------------------------------------------------------------------------
mat3d FEMat3dLocalElementMap::operator () (const FEMaterialPoint& mp)
{
	FEMesh& mesh = GetFEModel()->GetMesh();

	FEElement& el = *mp.m_elem;

	vec3d r0[FEElement::MAX_NODES];
	for (int i=0; i<el.Nodes(); ++i) r0[i] = mesh.Node(el.m_node[i]).m_r0;

	vec3d a, b, c, d;
	mat3d Q;

	a = r0[m_n[1] - 1] - r0[m_n[0] - 1];
	a.unit();

	if (m_n[2] != m_n[1])
	{
		d = r0[m_n[2] - 1] - r0[m_n[0] - 1];
	}
	else
	{
		d = vec3d(0,1,0);
		if (fabs(d*a) > 0.999) d = vec3d(1,0,0);
	}

	c = a^d;
	b = c^a;

	b.unit();
	c.unit();

	Q[0][0] = a.x; Q[0][1] = b.x; Q[0][2] = c.x;
	Q[1][0] = a.y; Q[1][1] = b.y; Q[1][2] = c.y;
	Q[2][0] = a.z; Q[2][1] = b.z; Q[2][2] = c.z;

	return Q;
}

//-----------------------------------------------------------------------------
FEMat3dValuator* FEMat3dLocalElementMap::copy()
{
	FEMat3dLocalElementMap* map = new FEMat3dLocalElementMap(GetFEModel());
	map->m_n[0] = m_n[0];
	map->m_n[1] = m_n[1];
	map->m_n[2] = m_n[2];
	return map;
}

//-----------------------------------------------------------------------------
void FEMat3dLocalElementMap::Serialize(DumpStream& ar)
{
	if (ar.IsSaving())
	{
		ar << m_n[0] << m_n[1] << m_n[2];
	}
	else
	{
		ar >> m_n[0] >> m_n[1] >> m_n[2];
	}
}

//=============================================================================
// FESphericalMap
//-----------------------------------------------------------------------------

BEGIN_FECORE_CLASS(FEMat3dSphericalMap, FEMat3dValuator)
	ADD_PARAMETER(m_c, "center");
	ADD_PARAMETER(m_r, "vector");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEMat3dSphericalMap::FEMat3dSphericalMap(FEModel* pfem): FEMat3dValuator(pfem)
{
	m_c = vec3d(0,0,0);
	m_r = vec3d(1,0,0);
}

//-----------------------------------------------------------------------------
bool FEMat3dSphericalMap::Init()
{
	return true;
}

//-----------------------------------------------------------------------------
mat3d FEMat3dSphericalMap::operator () (const FEMaterialPoint& mp)
{
	vec3d a = mp.m_r0;
	a -= m_c;
	a.unit();

	// setup the rotation vector
	vec3d x_unit(1,0,0);
	quatd q(x_unit, a);

	vec3d v = m_r;
	v.unit();
	q.RotateVector(v);
	a = v;

	vec3d d(0,1,0);
	d.unit();
	if (fabs(a*d) > .99) 
	{
		d = vec3d(0,0,1);
		d.unit();
	}

	vec3d c = a^d;
	vec3d b = c^a;

	a.unit();
	b.unit();
	c.unit();

	mat3d Q;
	Q[0][0] = a.x; Q[0][1] = b.x; Q[0][2] = c.x;
	Q[1][0] = a.y; Q[1][1] = b.y; Q[1][2] = c.y;
	Q[2][0] = a.z; Q[2][1] = b.z; Q[2][2] = c.z;

	return Q;
}

//-----------------------------------------------------------------------------
FEMat3dValuator* FEMat3dSphericalMap::copy()
{
	FEMat3dSphericalMap* map = new FEMat3dSphericalMap(GetFEModel());
	map->m_c = m_c;
	map->m_r = m_r;
	return map;
}

//=============================================================================
// FECylindricalMap
//-----------------------------------------------------------------------------

BEGIN_FECORE_CLASS(FEMat3dCylindricalMap, FEMat3dValuator)
	ADD_PARAMETER(m_c, "center");
	ADD_PARAMETER(m_a, "axis"  );
	ADD_PARAMETER(m_r, "vector");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEMat3dCylindricalMap::FEMat3dCylindricalMap(FEModel* pfem) : FEMat3dValuator(pfem)
{
	m_c = vec3d(0,0,0);
	m_a = vec3d(0,0,1);
	m_r = vec3d(1,0,0);
}

//-----------------------------------------------------------------------------
bool FEMat3dCylindricalMap::Init()
{
	m_a.unit();
	m_r.unit();
	return true;
}

//-----------------------------------------------------------------------------
mat3d FEMat3dCylindricalMap::operator () (const FEMaterialPoint& mp)
{
	// get the position of the material point
	vec3d p = mp.m_r0;

	// find the vector to the axis
	vec3d b = (p - m_c) - m_a*(m_a*(p - m_c)); b.unit();

	// setup the rotation vector
	vec3d x_unit(vec3d(1,0,0));
	quatd q(x_unit, b);

	// rotate the reference vector
	vec3d r(m_r); r.unit();
	q.RotateVector(r);

	// setup a local coordinate system with r as the x-axis
	vec3d d(0,1,0);
	q.RotateVector(d);
	if (fabs(d*r) > 0.99)
	{
		d = vec3d(0,0,1);
		q.RotateVector(d);
	}

	// find basis vectors
	vec3d e1 = r;
	vec3d e3 = (e1 ^ d); e3.unit();
	vec3d e2 = e3 ^ e1;

	// setup rotation matrix
	mat3d Q;
	Q[0][0] = e1.x; Q[0][1] = e2.x; Q[0][2] = e3.x;
	Q[1][0] = e1.y; Q[1][1] = e2.y; Q[1][2] = e3.y;
	Q[2][0] = e1.z; Q[2][1] = e2.z; Q[2][2] = e3.z;

	return Q;
}

//-----------------------------------------------------------------------------
FEMat3dValuator* FEMat3dCylindricalMap::copy()
{
	FEMat3dCylindricalMap* val = new FEMat3dCylindricalMap(GetFEModel());
	val->m_c = m_c;
	val->m_a = m_a;
	val->m_r = m_r;
	return val;
}

//=============================================================================
// FEPolarMap
//-----------------------------------------------------------------------------

BEGIN_FECORE_CLASS(FEMat3dPolarMap, FEMat3dValuator)
	ADD_PARAMETER(m_c, "center");
	ADD_PARAMETER(m_a, "axis"  );
	ADD_PARAMETER(m_d0, "vector1");
	ADD_PARAMETER(m_d1, "vector2");
	ADD_PARAMETER(m_R0, "radius1");
	ADD_PARAMETER(m_R1, "radius2");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEMat3dPolarMap::FEMat3dPolarMap(FEModel* pfem) : FEMat3dValuator(pfem)
{
	m_c = vec3d(0,0,0);
	m_a = vec3d(0,0,1);
	m_d0 = m_d1 = vec3d(1,0,0);
	m_R0 = 0; 
	m_R1 = 1;
}

//-----------------------------------------------------------------------------
bool FEMat3dPolarMap::Init()
{
	m_a.unit();
	m_d0.unit();
	m_d1.unit();
	return true;
}

//-----------------------------------------------------------------------------
mat3d FEMat3dPolarMap::operator () (const FEMaterialPoint& mp)
{
	// get the nodal position of material point
	vec3d p = mp.m_r0;

	// find the vector to the axis and its length
	vec3d b = (p - m_c) - m_a*(m_a*(p - m_c)); 
	double R = b.unit();

	// get the relative radius
	double R0 = m_R0;
	double R1 = m_R1;
	if (R1 == R0) R1 += 1;
	double w = (R - R0)/(R1 - R0);

	// get the fiber vectors
	vec3d v0 = m_d0;
	vec3d v1 = m_d1;
	quatd Q0(0,vec3d(0,0,1)), Q1(v0,v1);
	quatd Qw = quatd::slerp(Q0, Q1, w);
	vec3d v = v0; Qw.RotateVector(v);

	// setup the rotation vector
	vec3d x_unit(vec3d(1,0,0));
	quatd q(x_unit, b);

	// rotate the reference vector
	q.RotateVector(v);

	// setup a local coordinate system with r as the x-axis
	vec3d d(vec3d(0,1,0));
	q.RotateVector(d);
	if (fabs(d*v) > 0.99)
	{
		d = vec3d(0,0,1);
		q.RotateVector(d);
	}

	// find basis vectors
	vec3d e1 = v;
	vec3d e3 = (e1 ^ d); e3.unit();
	vec3d e2 = e3 ^ e1;

	// setup rotation matrix
	mat3d Q;
	Q[0][0] = e1.x; Q[0][1] = e2.x; Q[0][2] = e3.x;
	Q[1][0] = e1.y; Q[1][1] = e2.y; Q[1][2] = e3.y;
	Q[2][0] = e1.z; Q[2][1] = e2.z; Q[2][2] = e3.z;

	return Q;
}

//-----------------------------------------------------------------------------
FEMat3dValuator* FEMat3dPolarMap::copy()
{
	FEMat3dPolarMap* map = new FEMat3dPolarMap(GetFEModel());
	map->m_c = m_c;
	map->m_a = m_a;
	map->m_d0 = m_d0;
	map->m_d1 = m_d1;
	map->m_R0 = m_R0;
	map->m_R1 = m_R1;
	return map;
}

//=============================================================================
// FEVectorMap
//-----------------------------------------------------------------------------

BEGIN_FECORE_CLASS(FEMat3dVectorMap, FEMat3dValuator)
	ADD_PARAMETER(m_a, "a");
	ADD_PARAMETER(m_d, "d");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEMat3dVectorMap::FEMat3dVectorMap(FEModel* pfem) : FEMat3dValuator(pfem)
{
	m_a = vec3d(1,0,0);
	m_d = vec3d(0,1,0);
	m_Q.unit();
}

//-----------------------------------------------------------------------------
bool FEMat3dVectorMap::Init()
{
	// generators have to be unit vectors
	m_a.unit();
	m_d.unit();

	// make sure the vectors are not 0-vectors
	if ((m_a.norm2() == 0.0) || (m_d.norm2() == 0.0)) return false;

	// make sure that a, d are not aligned
	if (fabs(m_a*m_d) > 0.999)
	{
		// if they are, find a different value for d
		// Note: If this is used for a mat_axis parameter, then this
		// would modify the user specified value. 
		m_d = vec3d(1,0,0);
		if (fabs(m_a*m_d) > 0.999) m_d = vec3d(0,1,0);
	}

	vec3d a = m_a;
	vec3d d = m_d;

	vec3d c = a^d;
	vec3d b = c^a;

	a.unit();
	b.unit();
	c.unit();

	m_Q[0][0] = a.x; m_Q[0][1] = b.x; m_Q[0][2] = c.x;
	m_Q[1][0] = a.y; m_Q[1][1] = b.y; m_Q[1][2] = c.y;
	m_Q[2][0] = a.z; m_Q[2][1] = b.z; m_Q[2][2] = c.z;

	return true;
}

//-----------------------------------------------------------------------------
void FEMat3dVectorMap::SetVectors(vec3d a, vec3d d)
{ 
	m_a = a; 
	m_d = d; 
}

//-----------------------------------------------------------------------------
mat3d FEMat3dVectorMap::operator () (const FEMaterialPoint& mp)
{
	return m_Q;
}

//-----------------------------------------------------------------------------
FEMat3dValuator* FEMat3dVectorMap::copy()
{
	FEMat3dVectorMap* map = new FEMat3dVectorMap(GetFEModel());
	map->m_a = m_a;
	map->m_d = m_d;
	map->m_Q = m_Q;
	return map;
}

//-----------------------------------------------------------------------------
void FEMat3dVectorMap::Serialize(DumpStream &ar)
{
	if (ar.IsSaving())
	{
		ar << m_a << m_d << m_Q;
	}
	else
	{
		ar >> m_a >> m_d >> m_Q;
	}
}

//=============================================================================
// FESphericalAngleMap
//-----------------------------------------------------------------------------

BEGIN_FECORE_CLASS(FEMat3dSphericalAngleMap, FEMat3dValuator)
	ADD_PARAMETER(m_theta, "theta");
	ADD_PARAMETER(m_phi  , "phi"  );
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEMat3dSphericalAngleMap::FEMat3dSphericalAngleMap(FEModel* pfem) : FEMat3dValuator(pfem)
{
	m_theta = 0.0;
	m_phi   = 90.0;
	m_Q.unit();
}

//-----------------------------------------------------------------------------
bool FEMat3dSphericalAngleMap::Init()
{
	// convert from degress to radians
	const double pi = 4 * atan(1.0);
	const double the = m_theta*pi / 180.;
	const double phi = m_phi*pi / 180.;

	// define the first axis (i.e. the fiber vector)
	vec3d a;
	a.x = cos(the)*sin(phi);
	a.y = sin(the)*sin(phi);
	a.z = cos(phi);

	// define the second axis
	// and make sure it is not colinear with the first
	vec3d d(0, 0, 1);
	if (fabs(a*d) > 0.9) d = vec3d(0, 1, 0);

	// calculate the orthonormal axes
	vec3d c = a^d;
	vec3d b = c^a;
	a.unit();
	b.unit();
	c.unit();

	// setup the rotation matrix
	m_Q[0][0] = a.x; m_Q[0][1] = b.x; m_Q[0][2] = c.x;
	m_Q[1][0] = a.y; m_Q[1][1] = b.y; m_Q[1][2] = c.y;
	m_Q[2][0] = a.z; m_Q[2][1] = b.z; m_Q[2][2] = c.z;

	return true;
}

//-----------------------------------------------------------------------------
void FEMat3dSphericalAngleMap::SetAngles(double theta, double phi)
{ 
	m_theta = theta; 
	m_phi = phi; 
}

//-----------------------------------------------------------------------------
mat3d FEMat3dSphericalAngleMap::operator () (const FEMaterialPoint& mp)
{
	return m_Q;
}

//-----------------------------------------------------------------------------
FEMat3dValuator* FEMat3dSphericalAngleMap::copy()
{
	FEMat3dSphericalAngleMap* map = new FEMat3dSphericalAngleMap(GetFEModel());
	map->m_theta = m_theta;
	map->m_phi = m_phi;
	map->m_Q = m_Q;
	return map;
}

//-----------------------------------------------------------------------------
void FEMat3dSphericalAngleMap::Serialize(DumpStream &ar)
{
	if (ar.IsSaving())
	{
		ar << m_theta << m_phi << m_Q;
	}
	else
	{
		ar >> m_theta >> m_phi >> m_Q;
	}
}

//=============================================================================
// FEMappedValueMat3d
//-----------------------------------------------------------------------------

FEMappedValueMat3d::FEMappedValueMat3d(FEModel* fem) : FEMat3dValuator(fem)
{
	m_val = nullptr;
	m_scale = 0.0;
}

void FEMappedValueMat3d::setDataMap(FEDataMap* val, double scl)
{
	m_val = val;
	m_scale = scl;
}

mat3d FEMappedValueMat3d::operator()(const FEMaterialPoint& pt)
{
	return m_val->valueMat3d(pt)*m_scale;
}

FEMat3dValuator* FEMappedValueMat3d::copy()
{
	FEMappedValueMat3d* map = new FEMappedValueMat3d(GetFEModel());
	map->m_val = m_val;
	map->m_scale = m_scale;
	return map;
}
