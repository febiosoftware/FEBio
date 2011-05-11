#include "stdafx.h"
#include "FEPointBodyForce.h"
#include "fem.h"

BEGIN_PARAMETER_LIST(FEPointBodyForce, FEBodyForce);
	ADD_PARAMETER(m_a, FE_PARAM_DOUBLE, "a");
	ADD_PARAMETER(m_b, FE_PARAM_DOUBLE, "b");
	ADD_PARAMETER(m_rc, FE_PARAM_VEC3D, "rc");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEPointBodyForce::FEPointBodyForce(FEModel* pfem) : FEBodyForce(pfem)
{
	s[0] = s[1] = s[2] = 1.0; 
	m_rlc[0] = m_rlc[1] = m_rlc[2] = -1; 
	m_pel = 0; 
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
