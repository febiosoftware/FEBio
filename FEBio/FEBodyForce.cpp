#include "stdafx.h"
#include "FEBodyForce.h"

//-----------------------------------------------------------------------------
FENonConstBodyForce::FENonConstBodyForce()
{
	m_sz[0][0] = 0;
	m_sz[1][0] = 0;
	m_sz[2][0] = 0;
}

//-----------------------------------------------------------------------------
vec3d FENonConstBodyForce::force(FEMaterialPoint &mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// get the material point's spatial position
	vec3d r = pt.rt;

	// define a math parser object
	MathParser m;
	m.SetVariable("x", r.x);
	m.SetVariable("y", r.y);
	m.SetVariable("z", r.z);

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
