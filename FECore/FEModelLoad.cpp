#include "stdafx.h"
#include "FEModelLoad.h"
#include "FESolver.h"
#include "FERigidBody.h"

//-----------------------------------------------------------------------------
FEModelLoad::FEModelLoad(SUPER_CLASS_ID sid, FEModel* pfem) : FEModelComponent(sid, pfem)
{
}

//-----------------------------------------------------------------------------
void FEModelLoad::Residual(FEGlobalVector& R, const FETimePoint& pt)
{
}

//-----------------------------------------------------------------------------
void FEModelLoad::StiffnessMatrix(FESolver* psolver, const FETimePoint& pt)
{
}


//=============================================================================

BEGIN_PARAMETER_LIST(FERigidAxialForce, FEModelLoad);
	ADD_PARAMETER(m_ida, FE_PARAM_INT, "rbA");
	ADD_PARAMETER(m_idb, FE_PARAM_INT, "rbB");
	ADD_PARAMETER(m_ra0, FE_PARAM_VEC3D, "ra");
	ADD_PARAMETER(m_rb0, FE_PARAM_VEC3D, "rb");
	ADD_PARAMETER(m_s  , FE_PARAM_DOUBLE, "force");
	ADD_PARAMETER(m_brelative, FE_PARAM_BOOL, "relative");
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
	FEModel& fem = *GetFEModel();
	FERigidBody& bodyA = static_cast<FERigidBody&>(*fem.Object(m_ida));
	FERigidBody& bodyB = static_cast<FERigidBody&>(*fem.Object(m_idb));

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
	FEModel& fem = *GetFEModel();
	FERigidBody& bodyA = static_cast<FERigidBody&>(*fem.Object(m_ida));
	FERigidBody& bodyB = static_cast<FERigidBody&>(*fem.Object(m_idb));

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
	FEModel& fem = *GetFEModel();
	FERigidBody& bodyA = static_cast<FERigidBody&>(*fem.Object(m_ida));
	FERigidBody& bodyB = static_cast<FERigidBody&>(*fem.Object(m_idb));

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

	// build the stiffness matrix
	double M[3][3];
	M[0][0] = 1.0 - N.x*N.x; M[0][1] =     - N.x*N.y; M[0][2] =     - N.x*N.z;
	M[1][0] =     - N.y*N.x; M[1][1] = 1.0 - N.y*N.y; M[1][2] =     - N.y*N.z;
	M[2][0] =     - N.z*N.x; M[2][1] =     - N.z*N.y; M[2][2] = 1.0 - N.z*N.z;

	double A[3][3];
	A[0][0] =  0.0; A[0][1] = da.z; A[0][2] = -da.y;
	A[1][0] =-da.z; A[1][1] =  0.0; A[1][2] =  da.x;
	A[2][0] = da.y; A[2][1] =-da.x; A[2][2] =   0.0;

	double B[3][3];
	B[0][0] =  0.0; B[0][1] = db.z; B[0][2] = -db.y;
	B[1][0] =-db.z; B[1][1] =  0.0; B[1][2] =  db.x;
	B[2][0] = db.y; B[2][1] =-db.x; B[2][2] =   0.0;

	double S[3][3];
	S[0][0] = 0.0; S[0][1] = N.z; S[0][2] = -N.y;
	S[1][0] =-N.z; S[1][1] = 0.0; S[1][2] =  N.x;
	S[2][0] = N.y; S[2][1] =-N.x; S[2][2] =  0.0;

	double MA[3][3], MB[3][3];
	for (int i=0; i<3; ++i)
		for (int j=0; j<3; ++j)
		{
			MA[i][j] = 0.0;
			MB[i][j] = 0.0;
			for (int k=0; k<3; ++k)
			{
				MA[i][j] += M[i][k]*A[k][j];
				MB[i][j] += M[i][k]*B[k][j];
			}
		}

	double AMA[3][3], BMB[3][3], AMB[3][3];
	for (int i=0; i<3; ++i)
		for (int j=0; j<3; ++j)
		{
			AMA[i][j] = 0.0;
			BMB[i][j] = 0.0;
			AMB[i][j] = 0.0;
			for (int k=0; k<3; ++k)
			{
				AMA[i][j] += A[i][k]*MA[k][j];
				BMB[i][j] += B[i][k]*MB[k][j];
				AMB[i][j] += A[i][k]*MB[k][j];
			}
		}

	double SA[3][3], SB[3][3];
	for (int i=0; i<3; ++i)
		for (int j=0; j<3; ++j)
		{
			SA[i][j] = 0.0;
			SB[i][j] = 0.0;
			for (int k=0; k<3; ++k)
			{
				SA[i][j] += S[i][k]*A[k][j];
				SB[i][j] += S[i][k]*B[k][j];
			}
		}

	matrix K(12,12);
	K[0][0] = -M[0][0]; K[0][1] = -M[0][1]; K[0][2] = -M[0][2];
	K[1][0] = -M[1][0]; K[1][1] = -M[1][1]; K[1][2] = -M[1][2];
	K[2][0] = -M[2][0]; K[2][1] = -M[2][1]; K[2][2] = -M[2][2];

	K[0][3] = -MA[0][0]; K[0][4] = -MA[0][1]; K[0][5] = -MA[0][2];
	K[1][3] = -MA[1][0]; K[1][4] = -MA[1][1]; K[1][5] = -MA[1][2];
	K[2][3] = -MA[2][0]; K[2][4] = -MA[2][1]; K[2][5] = -MA[2][2];

	K[0][6] = M[0][0]; K[0][7] = M[0][1]; K[0][8] = M[0][2];
	K[1][6] = M[1][0]; K[1][7] = M[1][1]; K[1][8] = M[1][2];
	K[2][6] = M[2][0]; K[2][7] = M[2][1]; K[2][8] = M[2][2];

	K[0][9] = MB[0][0]; K[0][10] = MB[0][1]; K[0][11] = MB[0][2];
	K[1][9] = MB[1][0]; K[1][10] = MB[1][1]; K[1][11] = MB[1][2];
	K[2][9] = MB[2][0]; K[2][10] = MB[2][1]; K[2][11] = MB[2][2];

	K[3][3] = L*SA[0][0] + AMA[0][0]; K[3][4] = L*SA[0][1] + AMA[0][1]; K[3][5] = L*SA[0][2] + AMA[0][2];
	K[4][3] = L*SA[1][0] + AMA[1][0]; K[4][4] = L*SA[1][1] + AMA[1][1]; K[4][5] = L*SA[1][2] + AMA[1][2];
	K[5][3] = L*SA[2][0] + AMA[2][0]; K[5][4] = L*SA[2][1] + AMA[2][1]; K[5][5] = L*SA[2][2] + AMA[2][2];

	K[3][6] = MA[0][0]; K[3][7] = MA[1][0]; K[3][8] = MA[2][0];
	K[4][6] = MA[0][1]; K[4][7] = MA[1][1]; K[4][8] = MA[2][1];
	K[5][6] = MA[0][2]; K[5][7] = MA[1][2]; K[5][8] = MA[2][2];

	K[3][9] = -AMB[0][0]; K[3][10] = -AMB[0][1]; K[3][11] = -AMB[0][2];
	K[4][9] = -AMB[1][0]; K[4][10] = -AMB[1][1]; K[4][11] = -AMB[1][2];
	K[5][9] = -AMB[2][0]; K[5][10] = -AMB[2][1]; K[5][11] = -AMB[2][2];

	K[6][6] = -M[0][0]; K[6][7] = -M[0][1]; K[6][8] = -M[0][2];
	K[7][6] = -M[1][0]; K[7][7] = -M[1][1]; K[7][8] = -M[1][2];
	K[8][6] = -M[2][0]; K[8][7] = -M[2][1]; K[8][8] = -M[2][2];

	K[6][9] = -MB[0][0]; K[6][10] = -MB[0][1]; K[6][11] = -MB[0][2];
	K[7][9] = -MB[1][0]; K[7][10] = -MB[1][1]; K[7][11] = -MB[1][2];
	K[8][9] = -MB[2][0]; K[8][10] = -MB[2][1]; K[8][11] = -MB[2][2];

	K[ 9][9] = -L*SB[0][0] + BMB[0][0]; K[ 9][10] = -L*SB[0][1] + BMB[0][1]; K[ 9][11] = -L*SB[0][2] + BMB[0][2];
	K[10][9] = -L*SB[1][0] + BMB[1][0]; K[10][10] = -L*SB[1][1] + BMB[1][1]; K[10][11] = -L*SB[1][2] + BMB[1][2];
	K[11][9] = -L*SB[2][0] + BMB[2][0]; K[11][10] = -L*SB[2][1] + BMB[2][1]; K[11][11] = -L*SB[2][2] + BMB[2][2];

	// since this is a symmetric matrix, fill the bottom triangular part
	// and multiply by f
	for (int i=0; i<12; ++i)
	{
		K(i,i) *= f;
		for (int j=i+1; j<12; ++j)
		{
			K(i,j) *= f;
			K(j,i) = K(i,j);
		}
	}

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


//-----------------------------------------------------------------------------
FERigidBodyForce::FERigidBodyForce(FEModel* pfem) : FEModelLoad(FEBC_ID, pfem)
{
	m_bfollow = false;
}

//-----------------------------------------------------------------------------
void FERigidBodyForce::Serialize(DumpFile& ar)
{
	FEModelLoad::Serialize(ar);
	if (ar.IsSaving())
	{
		ar << ntype << bc << id << lc << sf;
	}
	else
	{
		ar >> ntype >> bc >> id >> lc >> sf;
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
	FERigidBody& rb = static_cast<FERigidBody&>(*fem.Object(id));

	if (m_bfollow == false)
	{
		int I  = rb.m_LM[bc];
		if (ntype == 0)
		{
			if ((I>=0) && (lc >= 0))
			{
				double f = fem.GetLoadCurve(lc)->Value()*sf;
				R[I] += f;
			}
		}
		else if (ntype == 1)
		{
			// get the current analysis step
			FEAnalysis* pstep = fem.GetCurrentStep();

			double t0 = pstep->m_tstart;
			double t1 = pstep->m_tend;
			double w = (tp.t - t0)/(t1 - t0);
			assert((w>=-0.0000001)&&(w<=1.0000001));
			double f0 = 0.0, f1 = sf;
			switch (bc)
			{
			case 0: f0 = rb.m_Fp.x; break;
			case 1: f0 = rb.m_Fp.y; break;
			case 2: f0 = rb.m_Fp.z; break;
			case 3: f0 = rb.m_Mp.x; break;
			case 4: f0 = rb.m_Mp.y; break;
			case 5: f0 = rb.m_Mp.z; break;
			}

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
