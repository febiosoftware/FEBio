#include "stdafx.h"
#include "FERigidForce.h"
#include "FECore/FERigidBody.h"
#include "FECore/FEModel.h"

//=============================================================================
BEGIN_PARAMETER_LIST(FERigidAxialForce, FEModelLoad);
	ADD_PARAMETER(m_ida      , FE_PARAM_INT   , "rbA"     );
	ADD_PARAMETER(m_idb      , FE_PARAM_INT   , "rbB"     );
	ADD_PARAMETER(m_ra0      , FE_PARAM_VEC3D , "ra"      );
	ADD_PARAMETER(m_rb0      , FE_PARAM_VEC3D , "rb"      );
	ADD_PARAMETER(m_s        , FE_PARAM_DOUBLE, "force"   );
	ADD_PARAMETER(m_brelative, FE_PARAM_BOOL  , "relative");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FERigidAxialForce::FERigidAxialForce(FEModel* pfem) : FEModelLoad(FEBC_ID, pfem)
{
	m_ida = m_idb = -1;
	m_ra0 = m_rb0 = vec3d(0,0,0);
	m_s = 0.0;
	m_brelative = false;
}

//-----------------------------------------------------------------------------
//! do some sanity checks
bool FERigidAxialForce::Init()
{
	// At this point the rigid ID's are still associated with the materials.
	// We want to associate them with the rigid objects instead.
	FEModel& fem = *GetFEModel();
	FEMaterial* pm = fem.GetMaterial(m_ida-1);
	m_ida = pm->GetRigidBodyID(); if (m_ida < 0) return false;
	pm = fem.GetMaterial(m_idb-1);
	m_idb = pm->GetRigidBodyID(); if (m_idb < 0) return false;

	// get the actual rigid bodies
	FERigidSystem& rigid = *GetFEModel()->GetRigidSystem();
	FERigidBody& bodyA = *rigid.Object(m_ida);
	FERigidBody& bodyB = *rigid.Object(m_idb);

	// get the attachment position in global coordinates for body A
	vec3d da0 = (m_brelative ? m_ra0 : m_ra0 - bodyA.m_r0);
	vec3d da = bodyA.m_qt*da0;
	vec3d a = da + bodyA.m_rt;

	// get the attachment position in global coordinates for body B
	vec3d db0 = (m_brelative ? m_rb0 : m_rb0 - bodyB.m_r0);
	vec3d db = bodyB.m_qt*db0;
	vec3d b = db + bodyB.m_rt;

	// get the unit axial vector
	// and make sure it is of finite length
	vec3d N = b - a; 
	double L = N.unit();
	if (L < 1e-17) return false;

	// all is well in the world
	return true;
}

//-----------------------------------------------------------------------------
void FERigidAxialForce::Serialize(DumpFile& ar)
{
	FEModelLoad::Serialize(ar);
	if (ar.IsSaving())
	{
		ar << m_ida << m_idb;
		ar << m_ra0 << m_rb0;
		ar << m_s << m_brelative;
	}
	else
	{
		ar >> m_ida >> m_idb;
		ar >> m_ra0 >> m_rb0;
		ar >> m_s >> m_brelative;
	}
}

//-----------------------------------------------------------------------------
//! Residual
void FERigidAxialForce::Residual(FEGlobalVector& R, const FETimePoint& tp)
{
	FERigidSystem& rigid = *GetFEModel()->GetRigidSystem();
	FERigidBody& bodyA = *rigid.Object(m_ida);
	FERigidBody& bodyB = *rigid.Object(m_idb);

	// get the attachment position in global coordinates for body A
	vec3d da0 = (m_brelative ? m_ra0 : m_ra0 - bodyA.m_r0);
	vec3d da = bodyA.m_qt*da0;
	vec3d a = da + bodyA.m_rt;

	// get the attachment position in global coordinates for body B
	vec3d db0 = (m_brelative ? m_rb0 : m_rb0 - bodyB.m_r0);
	vec3d db = bodyB.m_qt*db0;
	vec3d b = db + bodyB.m_rt;

	// get the unit axial vector
	vec3d N = b - a; N.unit();

	// calculate the force value
	double f = m_s;

	// get the axial force and torques
	vec3d F = N*f;
	vec3d Ma = da^F;
	vec3d Mb = db^F;

	// apply force and torque to body A
	int n;
	n = bodyA.m_LM[0]; if (n >= 0) R[n] += F.x;
	n = bodyA.m_LM[1]; if (n >= 0) R[n] += F.y;
	n = bodyA.m_LM[2]; if (n >= 0) R[n] += F.z;
	n = bodyA.m_LM[3]; if (n >= 0) R[n] += Ma.x;
	n = bodyA.m_LM[4]; if (n >= 0) R[n] += Ma.y;
	n = bodyA.m_LM[5]; if (n >= 0) R[n] += Ma.z;

	// apply force and torque to body B
	n = bodyB.m_LM[0]; if (n >= 0) R[n] -= F.x;
	n = bodyB.m_LM[1]; if (n >= 0) R[n] -= F.y;
	n = bodyB.m_LM[2]; if (n >= 0) R[n] -= F.z;
	n = bodyB.m_LM[3]; if (n >= 0) R[n] -= Mb.x;
	n = bodyB.m_LM[4]; if (n >= 0) R[n] -= Mb.y;
	n = bodyB.m_LM[5]; if (n >= 0) R[n] -= Mb.z;
}

//-----------------------------------------------------------------------------
//! Stiffness matrix
//! TODO: Only the stiffness contribution in the were the axial forces are applied
//!       to the center of mass has been implemented. 
void FERigidAxialForce::StiffnessMatrix(FESolver* psolver, const FETimePoint& tp)
{
	// Get the rigid bodies
	FERigidSystem& rigid = *GetFEModel()->GetRigidSystem();
	FERigidBody& bodyA = *rigid.Object(m_ida);
	FERigidBody& bodyB = *rigid.Object(m_idb);

	// get the attachment position in global coordinates for body A
	vec3d da0 = (m_brelative ? m_ra0 : m_ra0 - bodyA.m_r0);
	vec3d da = bodyA.m_qt*da0;
	vec3d pa = da + bodyA.m_rt;

	// get the attachment position in global coordinates for body B
	vec3d db0 = (m_brelative ? m_rb0 : m_rb0 - bodyB.m_r0);
	vec3d db = bodyB.m_qt*db0;
	vec3d pb = db + bodyB.m_rt;

	// setup the axial unit vector
	vec3d N = pb - pa; 
	double L = N.unit();

	// calculate the force value
	double f = - m_s / L;

	// build the stiffness matrix components
	mat3ds M = mat3dd(1.0) - dyad(N);
	mat3da A(-da);
	mat3da B(-db);
	mat3da S(-N);

	mat3d MA = M*A;
	mat3d MB = M*B;

	mat3d AMA = A*MA;
	mat3d BMB = B*MB;
	mat3d AMB = A*MB;

	mat3d SA = S*A;
	mat3d SB = S*B;

	// put it all together
	matrix K(12,12); K.zero();
	K.sub(0,0, M);
	K.sub(0,3, MA);
	K.add(0,6, M);
	K.add(0,9, MB);
	K.add(3,3, SA*L + AMA);
	K.add(3,6, MA);
	K.sub(3,9, AMB);
	K.sub(6,6,M);
	K.sub(6,9, MB);
	K.add(9,9, SB*(-L) + BMB);

	// since this is a symmetric matrix, fill the bottom triangular part
	K.copy_ut();

	// and multiply by f
	K *= f;

	// get the equation numbers
	vector<int> lm(12);
	lm[ 0] = bodyA.m_LM[0];
	lm[ 1] = bodyA.m_LM[1];
	lm[ 2] = bodyA.m_LM[2];
	lm[ 3] = bodyA.m_LM[3];
	lm[ 4] = bodyA.m_LM[4];
	lm[ 5] = bodyA.m_LM[5];

	lm[ 6] = bodyB.m_LM[0];
	lm[ 7] = bodyB.m_LM[1];
	lm[ 8] = bodyB.m_LM[2];
	lm[ 9] = bodyB.m_LM[3];
	lm[10] = bodyB.m_LM[4];
	lm[11] = bodyB.m_LM[5];

	// assemble into global matrix
	psolver->AssembleStiffness(lm, K);
}

//=============================================================================
FERigidBodyForce::FERigidBodyForce(FEModel* pfem) : FEModelLoad(FEBC_ID, pfem)
{
	m_bfollow = false;
	m_ntype = RAMP;
	m_trg = 0.0;
}

//-----------------------------------------------------------------------------
bool FERigidBodyForce::Init()
{
	// At this point the rigid ID's are still associated with the materials.
	// We want to associate them with the rigid objects instead.
	FEModel& fem = *GetFEModel();
	FEMaterial* pm = fem.GetMaterial(id-1);
	id = pm->GetRigidBodyID(); if (id < 0) return false;
	return true;
}

//-----------------------------------------------------------------------------
void FERigidBodyForce::Activate()
{
	FEModelLoad::Activate();

	FEModel& fem = *GetFEModel();
	FERigidSystem& rigid = *fem.GetRigidSystem();
	FERigidBody& rb = *rigid.Object(id);
	if (m_ntype == TARGET)
	{
		switch (bc)
		{
		case 0: m_trg = rb.m_Fr.x; break;
		case 1: m_trg = rb.m_Fr.y; break;
		case 2: m_trg= rb.m_Fr.z; break;
		case 3: m_trg = rb.m_Mr.x; break;
		case 4: m_trg = rb.m_Mr.y; break;
		case 5: m_trg= rb.m_Mr.z; break;
		}
	}
}

//-----------------------------------------------------------------------------
void FERigidBodyForce::Serialize(DumpFile& ar)
{
	FEModelLoad::Serialize(ar);
	if (ar.IsSaving())
	{
		ar << m_ntype << bc << id << lc << sf << m_trg;
	}
	else
	{
		ar >> m_ntype >> bc >> id >> lc >> sf >> m_trg;
	}
}

//-----------------------------------------------------------------------------
double FERigidBodyForce::Value()
{
	FEModel& fem = *GetFEModel();
	if (lc >= 0)
		return fem.GetLoadCurve(lc)->Value()*sf;
	else
		return sf;
}

//-----------------------------------------------------------------------------
//! Residual
void FERigidBodyForce::Residual(FEGlobalVector& R, const FETimePoint& tp)
{
	FEModel& fem = *GetFEModel();
	FERigidSystem& rigid = *fem.GetRigidSystem();
	FERigidBody& rb = *rigid.Object(id);

	if (m_bfollow == false)
	{
		int I  = rb.m_LM[bc];
		if (m_ntype == RAMP)
		{
			if ((I>=0) && (lc >= 0))
			{
				double f = fem.GetLoadCurve(lc)->Value()*sf;
				R[I] += f;
			}
		}
		else if (m_ntype == TARGET)
		{
			// get the current analysis step
			FEAnalysis* pstep = fem.GetCurrentStep();

			double t0 = pstep->m_tstart;
			double t1 = pstep->m_tend;
			double w = (tp.t - t0)/(t1 - t0);
			assert((w>=-0.0000001)&&(w<=1.0000001));
			double f0 = m_trg, f1 = sf;

			double f = f0*(1.0 - w) + f1*w;
			R[I] += f;
		}
	}
	else
	{
		// setup the force vector
		vec3d f(0,0,0);
		if      (bc == 0) f.x = Value();
		else if (bc == 1) f.y = Value();
		else if (bc == 2) f.z = Value();
	
		// apply the rigid body rotation
		f = rb.m_qt*f;

		// add to the residual
		int n;
		n = rb.m_LM[0]; if (n >= 0) R[n] += f.x;
		n = rb.m_LM[1]; if (n >= 0) R[n] += f.y;
		n = rb.m_LM[2]; if (n >= 0) R[n] += f.z;
	}
}

//-----------------------------------------------------------------------------
//! Stiffness matrix
void FERigidBodyForce::StiffnessMatrix(FESolver* psolver, const FETimePoint& tp)
{
	// I think for follower forces, I need to contribute to the stiffness matrix, but I'm not sure yet what.
}
