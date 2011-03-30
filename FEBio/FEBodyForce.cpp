#include "stdafx.h"
#include "FEBodyForce.h"

//=============================================================================
void FEBodyForce::Serialize(DumpFile& ar)
{
	if (ar.IsSaving())
	{
		ar << s[0] << s[1] << s[2];
		ar << lc[0] << lc[1] << lc[2];
	}
	else
	{
		ar >> s[0] >> s[1] >> s[2];
		ar >> lc[0] >> lc[1] >> lc[2];
	}
}

//=============================================================================

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

void FECentrifugalBodyForce::Serialize(DumpFile &ar)
{
	FEBodyForce::Serialize(ar);
	if (ar.IsSaving())
	{
		ar << c << n;
	}
	else
	{
		ar >> c >> n;
	}
}
