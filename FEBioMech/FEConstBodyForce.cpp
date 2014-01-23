#include "stdafx.h"
#include "FEConstBodyForce.h"

//=============================================================================
// FEConstBodyForce
//-----------------------------------------------------------------------------
BEGIN_PARAMETER_LIST(FEConstBodyForce, FEBodyForce);
	ADD_PARAMETER(m_f.x, FE_PARAM_DOUBLE, "x");
	ADD_PARAMETER(m_f.y, FE_PARAM_DOUBLE, "y");
	ADD_PARAMETER(m_f.z, FE_PARAM_DOUBLE, "z");
END_PARAMETER_LIST();

//=============================================================================
// FENonConstBodyForce
//-----------------------------------------------------------------------------
BEGIN_PARAMETER_LIST(FENonConstBodyForce, FEBodyForce);
	ADD_PARAMETER(m_sz[0], FE_PARAM_STRING, "x");
	ADD_PARAMETER(m_sz[1], FE_PARAM_STRING, "y");
	ADD_PARAMETER(m_sz[2], FE_PARAM_STRING, "z");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FENonConstBodyForce::FENonConstBodyForce(FEModel* pfem) : FEBodyForce(pfem)
{
	m_sz[0][0] = 0;
	m_sz[1][0] = 0;
	m_sz[2][0] = 0;
}

//-----------------------------------------------------------------------------
void FEConstBodyForce::Serialize(DumpFile &ar)
{
	FEBodyForce::Serialize(ar);
}

//-----------------------------------------------------------------------------
vec3d FENonConstBodyForce::force(FEMaterialPoint &mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// get the material point's spatial position
	vec3d r = pt.m_rt;

	// define a math parser object
	MathParser m;
	m.SetVariable("x", r.x);
	m.SetVariable("y", r.y);
	m.SetVariable("z", r.z);
	m.SetVariable("t", mp.time);

	// calculate the force
	vec3d f;
	int ierr;
	f.x = m.eval(m_sz[0], ierr);
	f.y = m.eval(m_sz[1], ierr);
	f.z = m.eval(m_sz[2], ierr);
	
	return f;
}

//-----------------------------------------------------------------------------
mat3ds FENonConstBodyForce::stiffness(FEMaterialPoint& pt)
{
	return mat3ds(0,0,0,0,0,0);
}

//-----------------------------------------------------------------------------
void FENonConstBodyForce::Serialize(DumpFile &ar)
{
	FEBodyForce::Serialize(ar);
	if (ar.IsSaving())
	{
		ar << m_sz[0] << m_sz[1] << m_sz[2];
	}
	else
	{
		ar >> m_sz[0] >> m_sz[1] >> m_sz[2];
	}
}

//=============================================================================
// FECentrifugalBodyForce
//-----------------------------------------------------------------------------
BEGIN_PARAMETER_LIST(FECentrifugalBodyForce, FEBodyForce);
	ADD_PARAMETER(w, FE_PARAM_DOUBLE, "angular_speed");
	ADD_PARAMETER(n, FE_PARAM_VEC3D, "rotation_axis");
	ADD_PARAMETER(c, FE_PARAM_VEC3D, "rotation_center");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
void FECentrifugalBodyForce::Serialize(DumpFile &ar)
{
	FEBodyForce::Serialize(ar);
	if (ar.IsSaving())
	{
		ar << w << c << n;
	}
	else
	{
		ar >> w >> c >> n;
	}
}
