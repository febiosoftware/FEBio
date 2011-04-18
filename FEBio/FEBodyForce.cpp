#include "stdafx.h"
#include "FEBodyForce.h"
#include "fem.h"

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

//=============================================================================
// FEPointBodyForce
//=============================================================================

//-----------------------------------------------------------------------------
FEPointBodyForce::FEPointBodyForce(FEM* pfem)
{
	s[0] = s[1] = s[2] = 1.0; 
	m_rlc[0] = m_rlc[1] = m_rlc[2] = -1; 
	m_pel = 0; 
	m_pfem = pfem; 
	m_brigid = true; 
	m_ntype = POINT; 
	m_inode = -1; 
}

//-----------------------------------------------------------------------------
vec3d FEPointBodyForce::force(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	vec3d x = pt.rt;
	vec3d n = x - m_rc;
	double l = n.unit();

	double g = m_a*exp(-m_b*l);
	return n*g;
}

//-----------------------------------------------------------------------------
mat3ds FEPointBodyForce::stiffness(FEMaterialPoint &mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	vec3d x = pt.rt;
	vec3d n = x - m_rc;
	double l = n.unit();

	mat3ds k;
	if (l == 0.0)
	{
		k.zero();
	}
	else
	{
		double g = m_a*exp(-m_b*l);

		mat3ds nxn = dyad(n);
		mat3ds I = mat3dd(1.0);

		k = (nxn*m_b - (I - nxn)/l)*g;
	}

	return k;
}

//-----------------------------------------------------------------------------
void FEPointBodyForce::Serialize(DumpFile &ar)
{
	if (ar.IsSaving())
	{
		ar << m_a << m_b << m_rc;
		ar << m_rlc[0] << m_rlc[1] << m_rlc[2];
		ar << m_ntype << m_inode << m_brigid;
	}
	else
	{
		ar >> m_a >> m_b >> m_rc;
		ar >> m_rlc[0] >> m_rlc[1] >> m_rlc[2];
		ar >> m_ntype >> m_inode >> m_brigid;
	}
}

//-----------------------------------------------------------------------------
void FEPointBodyForce::Init()
{
	assert(m_pfem);

	if (m_ntype == POINT)
	{
		if (!m_brigid)
		{
			// make sure we don't move the point
			m_rlc[0] = m_rlc[1] = m_rlc[2] = -1;

			// find the element in which point r0 lies
			FEMesh& m = m_pfem->m_mesh;
			m_pel = m.FindSolidElement(m_rc, m_rs);
		}
		else m_pel = 0;
	}
	else 
	{
		// make sure we don't move the point
		m_rlc[0] = m_rlc[1] = m_rlc[2] = -1;

		assert(m_inode >= 0);
		FEMesh& m = m_pfem->m_mesh;
		m_rc = m.Node(m_inode).m_r0;
	}
}

//-----------------------------------------------------------------------------
// Update the position of the body force
void FEPointBodyForce::Update()
{
	if (m_ntype == POINT)
	{
		if (m_pel)
		{
			FEMesh& m = m_pfem->m_mesh;
			vec3d x[8];
			for (int i=0; i<8; ++i) x[i] = m.Node(m_pel->m_node[i]).m_rt;

			double* r = m_rs;
			double H[8];
			H[0] = 0.125*(1 - r[0])*(1 - r[1])*(1 - r[2]);
			H[1] = 0.125*(1 + r[0])*(1 - r[1])*(1 - r[2]);
			H[2] = 0.125*(1 + r[0])*(1 + r[1])*(1 - r[2]);
			H[3] = 0.125*(1 - r[0])*(1 + r[1])*(1 - r[2]);
			H[4] = 0.125*(1 - r[0])*(1 - r[1])*(1 + r[2]);
			H[5] = 0.125*(1 + r[0])*(1 - r[1])*(1 + r[2]);
			H[6] = 0.125*(1 + r[0])*(1 + r[1])*(1 + r[2]);
			H[7] = 0.125*(1 - r[0])*(1 + r[1])*(1 + r[2]);

			m_rc = vec3d(0,0,0);
			for (int i=0; i<8; ++i) m_rc += x[i]*H[i];
		}
	}
	else
	{
		FEMesh& m = m_pfem->m_mesh;
		m_rc = m.Node(m_inode).m_rt;
	}
}
