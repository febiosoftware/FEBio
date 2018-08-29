// FECoordSysMap.cpp: implementation of the FECoordSysMap class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FECoordSysMap.h"
#include "FEMesh.h"
#include "FEModel.h"
#include "DumpStream.h"
#include "FEElement.h"

//-----------------------------------------------------------------------------
FECoordSysMap::FECoordSysMap(FEModel* pfem) : FECoreBase(FECOORDSYSMAP_ID) 
{ 
	m_pfem = pfem; 
}

//-----------------------------------------------------------------------------
FECoordSysMap::~FECoordSysMap() {}

//-----------------------------------------------------------------------------
//! initialization
bool FECoordSysMap::Init() { return true; }


//=============================================================================
// FELocalMap
//-----------------------------------------------------------------------------

BEGIN_PARAMETER_LIST(FELocalMap, FECoordSysMap)
	ADD_PARAMETERV(m_n, FE_PARAM_INT, 3, "local");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FELocalMap::FELocalMap(FEModel* pfem) : FECoordSysMap(pfem)
{
	m_n[0] = -1;
	m_n[1] = -1;
	m_n[2] = -1;
}

//-----------------------------------------------------------------------------
bool FELocalMap::Init()
{
	if ((m_n[0]==-1)&&(m_n[1]==-1)&&(m_n[2]==-1)) { m_n[0] = 0; m_n[1] = 1; m_n[2] = 3; }
	if (m_n[2] == -1) m_n[2] = m_n[1];
	return true;
}

//-----------------------------------------------------------------------------
void FELocalMap::SetLocalNodes(int n1, int n2, int n3)
{
	m_n[0] = n1;
	m_n[1] = n2;
	m_n[2] = n3;
}

//-----------------------------------------------------------------------------
mat3d FELocalMap::LocalElementCoord(FEElement& el, int n)
{
	FEMesh& mesh = GetFEModel()->GetMesh();
	vec3d r0[FEElement::MAX_NODES];
	for (int i=0; i<el.Nodes(); ++i) r0[i] = mesh.Node(el.m_node[i]).m_r0;

	vec3d a, b, c, d;
	mat3d Q;

	a = r0[m_n[1]] - r0[m_n[0]];
	a.unit();

	if (m_n[2] != m_n[1])
	{
		d = r0[m_n[2]] - r0[m_n[0]];
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
void FELocalMap::Serialize(DumpStream& ar)
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

BEGIN_PARAMETER_LIST(FESphericalMap, FECoordSysMap)
	ADD_PARAMETER(m_c, FE_PARAM_VEC3D, "center");
	ADD_PARAMETER(m_r, FE_PARAM_VEC3D, "vector");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FESphericalMap::FESphericalMap(FEModel* pfem): FECoordSysMap(pfem)
{
	m_c = vec3d(0,0,0);
	m_r = vec3d(1,0,0);
}

//-----------------------------------------------------------------------------
bool FESphericalMap::Init()
{
	return true;
}

//-----------------------------------------------------------------------------
mat3d FESphericalMap::LocalElementCoord(FEElement& el, int n)
{
	FEMesh& mesh = GetFEModel()->GetMesh();
	vec3d r0[FEElement::MAX_NODES];
	for (int i=0; i<el.Nodes(); ++i) r0[i] = mesh.Node(el.m_node[i]).m_r0;

	double* H = el.H(n);
	vec3d a;
	for (int i=0; i<el.Nodes(); ++i) a += r0[i]*H[i];
	a -= m_c;
	a.unit();

	// setup the rotation vector
	vec3d x_unit(1,0,0);
	quatd q(x_unit, a);

	vec3d v = m_r;
	v.unit();
	q.RotateVector(v);
	a = v;

	vec3d d = r0[1] - r0[0];
	d.unit();
	if (fabs(a*d) > .99) 
	{
		d = r0[2] - r0[1];
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
void FESphericalMap::Serialize(DumpStream& ar)
{
	if (ar.IsSaving())
	{
		ar << m_c;
	}
	else
	{
		ar >> m_c;
	}
}


//=============================================================================
// FECylindricalMap
//-----------------------------------------------------------------------------

BEGIN_PARAMETER_LIST(FECylindricalMap, FECoordSysMap)
	ADD_PARAMETER(m_c, FE_PARAM_VEC3D, "center");
	ADD_PARAMETER(m_a, FE_PARAM_VEC3D, "axis"  );
	ADD_PARAMETER(m_r, FE_PARAM_VEC3D, "vector");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FECylindricalMap::FECylindricalMap(FEModel* pfem) : FECoordSysMap(pfem)
{
	m_c = vec3d(0,0,0);
	m_a = vec3d(0,0,1);
	m_r = vec3d(1,0,0);
}

//-----------------------------------------------------------------------------
bool FECylindricalMap::Init()
{
	m_a.unit();
	m_r.unit();
	return true;
}

//-----------------------------------------------------------------------------
mat3d FECylindricalMap::LocalElementCoord(FEElement& el, int n)
{
	FEMesh& mesh = GetFEModel()->GetMesh();
	// get the element nodes
	vec3d r0[FEElement::MAX_NODES];
	for (int i=0; i<el.Nodes(); ++i) r0[i] = mesh.Node(el.m_node[i]).m_r0;

	// find the nodal position of the integration point n
	vec3d p = el.Evaluate(r0, n);

	// find the vector to the axis
	vec3d b = (p - m_c) - m_a*(m_a*(p - m_c)); b.unit();

	// setup the rotation vector
	vec3d x_unit(vec3d(1,0,0));
	quatd q(x_unit, b);

	// rotate the reference vector
	vec3d r(m_r); r.unit();
	q.RotateVector(r);

	// setup a local coordinate system with r as the x-axis
	vec3d d(vec3d(0,1,0));
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
void FECylindricalMap::Serialize(DumpStream& ar)
{
	if (ar.IsSaving())
	{
		ar << m_c << m_a << m_r;
	}
	else
	{
		ar >> m_c >> m_a >> m_r;
	}
}

//=============================================================================
// FEPolarMap
//-----------------------------------------------------------------------------

BEGIN_PARAMETER_LIST(FEPolarMap, FECoordSysMap)
	ADD_PARAMETER(m_c, FE_PARAM_VEC3D, "center");
	ADD_PARAMETER(m_a, FE_PARAM_VEC3D, "axis"  );
	ADD_PARAMETER(m_d0, FE_PARAM_VEC3D, "vector1");
	ADD_PARAMETER(m_d1, FE_PARAM_VEC3D, "vector2");
	ADD_PARAMETER(m_R0, FE_PARAM_DOUBLE, "radius1");
	ADD_PARAMETER(m_R1, FE_PARAM_DOUBLE, "radius2");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEPolarMap::FEPolarMap(FEModel* pfem) : FECoordSysMap(pfem)
{
	m_c = vec3d(0,0,0);
	m_a = vec3d(0,0,1);
	m_d0 = m_d1 = vec3d(1,0,0);
	m_R0 = 0; 
	m_R1 = 1;
}

//-----------------------------------------------------------------------------
bool FEPolarMap::Init()
{
	m_a.unit();
	m_d0.unit();
	m_d1.unit();
	return true;
}

//-----------------------------------------------------------------------------
mat3d FEPolarMap::LocalElementCoord(FEElement& el, int n)
{
	FEMesh& mesh = GetFEModel()->GetMesh();

	// get the element nodes
	vec3d r0[FEElement::MAX_NODES];
	for (int i=0; i<el.Nodes(); ++i) r0[i] = mesh.Node(el.m_node[i]).m_r0;

	// find the nodal position of the integration point n
	vec3d p = el.Evaluate(r0, n);

	// find the vector to the axis and its lenght
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
void FEPolarMap::Serialize(DumpStream& ar)
{
	if (ar.IsSaving())
	{
		ar << m_c << m_a << m_d0 << m_d1 << m_R0 << m_R1;
	}
	else
	{
		ar >> m_c >> m_a >> m_d0 >> m_d1 >> m_R0 >> m_R1;
	}
}

//=============================================================================
// FEVectorMap
//-----------------------------------------------------------------------------

BEGIN_PARAMETER_LIST(FEVectorMap, FECoordSysMap)
	ADD_PARAMETER(m_a, FE_PARAM_VEC3D, "vector");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEVectorMap::FEVectorMap(FEModel* pfem) : FECoordSysMap(pfem) 
{
	m_a = vec3d(1,0,0);
	m_d = vec3d(0,1,0);
}

//-----------------------------------------------------------------------------
bool FEVectorMap::Init()
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
	return true;
}

//-----------------------------------------------------------------------------
void FEVectorMap::SetVectors(vec3d a, vec3d d)
{ 
	m_a = a; 
	m_d = d; 
}

//-----------------------------------------------------------------------------
mat3d FEVectorMap::LocalElementCoord(FEElement& el, int n)
{
	vec3d a = m_a;
	vec3d d = m_d;

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
void FEVectorMap::Serialize(DumpStream &ar)
{
	if (ar.IsSaving())
	{
		ar << m_a << m_d;
	}
	else
	{
		ar >> m_a >> m_d;
	}
}

//=============================================================================
// FESphericalAngleMap
//-----------------------------------------------------------------------------

BEGIN_PARAMETER_LIST(FESphericalAngleMap, FECoordSysMap)
	ADD_PARAMETER(m_theta, FE_PARAM_DOUBLE, "theta");
	ADD_PARAMETER(m_phi  , FE_PARAM_DOUBLE, "phi"  );
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FESphericalAngleMap::FESphericalAngleMap(FEModel* pfem) : FECoordSysMap(pfem)
{
	m_theta = 0.0;
	m_phi   = 90.0;
}

//-----------------------------------------------------------------------------
bool FESphericalAngleMap::Init()
{
	return true;
}

//-----------------------------------------------------------------------------
void FESphericalAngleMap::SetAngles(double theta, double phi)
{ 
	m_theta = theta; 
	m_phi = phi; 
}

//-----------------------------------------------------------------------------
mat3d FESphericalAngleMap::LocalElementCoord(FEElement& el, int n)
{
	// convert from degress to radians
	const double pi = 4*atan(1.0);
	const double the = m_theta*pi/180.;
	const double phi = m_phi*pi/180.;

	// define the first axis (i.e. the fiber vector)
	vec3d a;
	a.x = cos(the)*sin(phi);
	a.y = sin(the)*sin(phi);
	a.z = cos(phi);

	// define the second axis
	// and make sure it is not colinear with the first
	vec3d d(0,0,1);
	if (fabs(a*d) > 0.9) d = vec3d(0,1,0);

	// calculate the orthonormal axes
	vec3d c = a^d;
	vec3d b = c^a;
	a.unit();
	b.unit();
	c.unit();

	// setup the rotation matrix
	mat3d Q;
	Q[0][0] = a.x; Q[0][1] = b.x; Q[0][2] = c.x;
	Q[1][0] = a.y; Q[1][1] = b.y; Q[1][2] = c.y;
	Q[2][0] = a.z; Q[2][1] = b.z; Q[2][2] = c.z;

	return Q;
}

//-----------------------------------------------------------------------------
void FESphericalAngleMap::Serialize(DumpStream &ar)
{
	if (ar.IsSaving())
	{
		ar << m_theta << m_phi;
	}
	else
	{
		ar >> m_theta >> m_phi;
	}
}
