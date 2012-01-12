// FERigidJoint.cpp: implementation of the FERigidJoint class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FERigidJoint.h"
#include "fem.h"
#include "FESolver.h"
#include "FEBioLib/log.h"
#include "FEAnalysisStep.h"

//-----------------------------------------------------------------------------
BEGIN_PARAMETER_LIST(FERigidJoint, FEParamContainer);
	ADD_PARAMETER(m_atol, FE_PARAM_DOUBLE, "tolerance");
	ADD_PARAMETER(m_eps , FE_PARAM_DOUBLE, "penalty"  );
	ADD_PARAMETER(m_nRBa, FE_PARAM_INT   , "body_a"   );
	ADD_PARAMETER(m_nRBb, FE_PARAM_INT   , "body_b"   );
	ADD_PARAMETER(m_q0  , FE_PARAM_VEC3D , "joint"    );
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FERigidJoint::FERigidJoint(FEM* pfem)
{
	static int count = 1;
	m_pfem = pfem;
	m_nID = count++;
}

//-----------------------------------------------------------------------------
FERigidJoint::~FERigidJoint()
{
	m_pfem = 0;
}

//-----------------------------------------------------------------------------
void FERigidJoint::JointForces(vector<double>& R)
{
	int i;
	double fe[6];

	FERigidBody& RBa = m_pfem->m_RB[m_nRBa];
	FERigidBody& RBb = m_pfem->m_RB[m_nRBb];

	vec3d a = m_qa0;
	RBa.m_qt.RotateVector(a);

	fe[0] = m_F.x;
	fe[1] = m_F.y;
	fe[2] = m_F.z;

	fe[3] = a.y*m_F.z-a.z*m_F.y;
	fe[4] = a.z*m_F.x-a.x*m_F.z;
	fe[5] = a.x*m_F.y-a.y*m_F.x;

	for (i=0; i<6; ++i) if (RBa.m_LM[i] >= 0) R[RBa.m_LM[i]] += fe[i];

	a = m_qb0;
	RBb.m_qt.RotateVector(a);

	fe[0] = -m_F.x;
	fe[1] = -m_F.y;
	fe[2] = -m_F.z;

	fe[3] = -a.y*m_F.z+a.z*m_F.y;
	fe[4] = -a.z*m_F.x+a.x*m_F.z;
	fe[5] = -a.x*m_F.y+a.y*m_F.x;

	for (i=0; i<6; ++i) if (RBb.m_LM[i] >= 0) R[RBb.m_LM[i]] += fe[i];
}

//-----------------------------------------------------------------------------
void FERigidJoint::JointStiffness()
{
	int j, k;

	vector<int> LM(12);
	matrix ke(12,12);
	ke.zero();
	vec3d a;

	double y1[3][3], y2[3][3], y11[3][3], y12[3][3], y22[3][3];

	FERigidBody& RBa = m_pfem->m_RB[m_nRBa];
	FERigidBody& RBb = m_pfem->m_RB[m_nRBb];

	a = m_qa0;
	RBa.m_qt.RotateVector(a);

	y1[0][0] =    0; y1[0][1] =  a.z; y1[0][2] = -a.y;
	y1[1][0] = -a.z; y1[1][1] =    0; y1[1][2] =  a.x;
	y1[2][0] =  a.y; y1[2][1] = -a.x; y1[2][2] =    0;

	a = m_qb0;
	RBb.m_qt.RotateVector(a);

	y2[0][0] =    0; y2[0][1] =  a.z; y2[0][2] = -a.y;
	y2[1][0] = -a.z; y2[1][1] =    0; y2[1][2] =  a.x;
	y2[2][0] =  a.y; y2[2][1] = -a.x; y2[2][2] =    0;

	for (j=0; j<3; ++j)
	{
		y11[j][0] = y1[0][j]*y1[0][0]+y1[1][j]*y1[1][0]+y1[2][j]*y1[2][0];
		y11[j][1] = y1[0][j]*y1[0][1]+y1[1][j]*y1[1][1]+y1[2][j]*y1[2][1];
		y11[j][2] = y1[0][j]*y1[0][2]+y1[1][j]*y1[1][2]+y1[2][j]*y1[2][2];

		y12[j][0] = y1[0][j]*y2[0][0]+y1[1][j]*y2[1][0]+y1[2][j]*y2[2][0];
		y12[j][1] = y1[0][j]*y2[0][1]+y1[1][j]*y2[1][1]+y1[2][j]*y2[2][1];
		y12[j][2] = y1[0][j]*y2[0][2]+y1[1][j]*y2[1][2]+y1[2][j]*y2[2][2];

		y22[j][0] = y2[0][j]*y2[0][0]+y2[1][j]*y2[1][0]+y2[2][j]*y2[2][0];
		y22[j][1] = y2[0][j]*y2[0][1]+y2[1][j]*y2[1][1]+y2[2][j]*y2[2][1];
		y22[j][2] = y2[0][j]*y2[0][2]+y2[1][j]*y2[1][2]+y2[2][j]*y2[2][2];
	}

	ke[0][0] = ke[1][1] = ke[2][2] =  1;
	ke[0][6] = ke[1][7] = ke[2][8] = -1;
	ke[6][0] = ke[7][1] = ke[8][2] = -1;
	ke[6][6] = ke[7][7] = ke[8][8] =  1;

	ke[0][3] = y1[0][0]; ke[0][4] = y1[0][1]; ke[0][5] = y1[0][2];
	ke[1][3] = y1[1][0]; ke[1][4] = y1[1][1]; ke[1][5] = y1[1][2];
	ke[2][3] = y1[2][0]; ke[2][4] = y1[2][1]; ke[2][5] = y1[2][2];

	ke[0][9] = -y2[0][0]; ke[0][10] = -y2[0][1]; ke[0][11] = -y2[0][2];
	ke[1][9] = -y2[1][0]; ke[1][10] = -y2[1][1]; ke[1][11] = -y2[1][2];
	ke[2][9] = -y2[2][0]; ke[2][10] = -y2[2][1]; ke[2][11] = -y2[2][2];

	ke[3][0] = y1[0][0]; ke[3][1] = y1[1][0]; ke[3][2] = y1[2][0];
	ke[4][0] = y1[0][1]; ke[4][1] = y1[1][1]; ke[4][2] = y1[2][1];
	ke[5][0] = y1[0][2]; ke[5][1] = y1[1][2]; ke[5][2] = y1[2][2];

	ke[3][3] = y11[0][0]; ke[3][4] = y11[0][1]; ke[3][5] = y11[0][2];
	ke[4][3] = y11[1][0]; ke[4][4] = y11[1][1]; ke[4][5] = y11[1][2];
	ke[5][3] = y11[2][0]; ke[5][4] = y11[2][1]; ke[5][5] = y11[2][2];

	ke[3][6] = -y1[0][0]; ke[3][7] = -y1[1][0]; ke[3][8] = -y1[2][0];
	ke[4][6] = -y1[0][1]; ke[4][7] = -y1[1][1]; ke[4][8] = -y1[2][1];
	ke[5][6] = -y1[0][2]; ke[5][7] = -y1[1][2]; ke[5][8] = -y1[2][2];

	ke[3][9] = -y12[0][0]; ke[3][10] = -y12[0][1]; ke[3][11] = -y12[0][2];
	ke[4][9] = -y12[1][0]; ke[4][10] = -y12[1][1]; ke[4][11] = -y12[1][2];
	ke[5][9] = -y12[2][0]; ke[5][10] = -y12[2][1]; ke[5][11] = -y12[2][2];

	ke[6][3] = -y1[0][0]; ke[6][4] = -y1[0][1]; ke[6][5] = -y1[0][2];
	ke[7][3] = -y1[1][0]; ke[7][4] = -y1[1][1]; ke[7][5] = -y1[1][2];
	ke[8][3] = -y1[2][0]; ke[8][4] = -y1[2][1]; ke[8][5] = -y1[2][2];

	ke[6][9] = y2[0][0]; ke[6][10] = y2[0][1]; ke[6][11] = y2[0][2];
	ke[7][9] = y2[1][0]; ke[7][10] = y2[1][1]; ke[7][11] = y2[1][2];
	ke[8][9] = y2[2][0]; ke[8][10] = y2[2][1]; ke[8][11] = y2[2][2];

	ke[ 9][0] = -y2[0][0]; ke[ 9][1] = -y2[1][0]; ke[ 9][2] = -y2[2][0];
	ke[10][0] = -y2[0][1]; ke[10][1] = -y2[1][1]; ke[10][2] = -y2[2][1];
	ke[11][0] = -y2[0][2]; ke[11][1] = -y2[1][2]; ke[11][2] = -y2[2][2];

	ke[ 9][3] = -y12[0][0]; ke[ 9][4] = -y12[1][0]; ke[ 9][5] = -y12[2][0];
	ke[10][3] = -y12[0][1]; ke[10][4] = -y12[1][1]; ke[10][5] = -y12[2][1];
	ke[11][3] = -y12[0][2]; ke[11][4] = -y12[1][2]; ke[11][5] = -y12[2][2];

	ke[ 9][6] = y2[0][0]; ke[ 9][7] = y2[1][0]; ke[ 9][8] = y2[2][0];
	ke[10][6] = y2[0][1]; ke[10][7] = y2[1][1]; ke[10][8] = y2[2][1];
	ke[11][6] = y2[0][2]; ke[11][7] = y2[1][2]; ke[11][8] = y2[2][2];

	ke[ 9][9] = y22[0][0]; ke[ 9][10] = y22[0][1]; ke[ 9][11] = y22[0][2];
	ke[10][9] = y22[1][0]; ke[10][10] = y22[1][1]; ke[10][11] = y22[1][2];
	ke[11][9] = y22[2][0]; ke[11][10] = y22[2][1]; ke[11][11] = y22[2][2];

	for (j=0; j<12; ++j)
		for (k=0; k<12; ++k)
		{
			ke[j][k] *= -m_eps;
		}

	for (j=0; j<6; ++j)
	{
		LM[j  ] = RBa.m_LM[j];
		LM[j+6] = RBb.m_LM[j];
	}

	FEAnalysisStep* pstep = dynamic_cast<FEAnalysisStep*>(m_pfem->m_pStep);
	pstep->m_psolver->m_pK->Assemble(ke, LM);
}

//-----------------------------------------------------------------------------
bool FERigidJoint::Augment()
{
	vec3d ra, rb, qa, qb, c,  Lm;
	double normF0, normF1;
	bool bconv = true;

	FERigidBody& RBa = m_pfem->m_RB[m_nRBa];
	FERigidBody& RBb = m_pfem->m_RB[m_nRBb];

	ra = RBa.m_rt;
	rb = RBb.m_rt;

	qa = m_qa0;
	RBa.m_qt.RotateVector(qa);

	qb = m_qb0;
	RBb.m_qt.RotateVector(qb);

	c = ra + qa - rb - qb;

	normF0 = sqrt(m_L*m_L);

	// calculate trial multiplier
	Lm = m_L + c*m_eps;

	normF1 = sqrt(Lm*Lm);

	// check convergence of constraints
	clog.printf(" rigid joint # %d\n", m_nID);
	clog.printf("                  CURRENT        REQUIRED\n");
	double pctn = 0;
	if (fabs(normF1) > 1e-10) pctn = fabs((normF1 - normF0)/normF1);
	clog.printf("    force : %15le %15le\n", pctn, m_atol);
	clog.printf("    gap   : %15le       ***\n", c.norm());
		
	if (pctn >= m_atol)
	{
		bconv = false;

		// update multiplier
		m_L = m_L + c*m_eps;
	
		// update force
		m_F = m_L + c*m_eps;
	}
	
	return bconv;
}

//-----------------------------------------------------------------------------
void FERigidJoint::Serialize(DumpFile& ar)
{
	if (ar.IsSaving())
	{
		ar << m_nID;
		ar << m_nRBa << m_nRBb;
		ar << m_q0 << m_qa0 << m_qb0;
		ar << m_F << m_L << m_eps << m_atol;
	}
	else
	{
		ar >> m_nID;
		ar >> m_nRBa >> m_nRBb;
		ar >> m_q0 >> m_qa0 >> m_qb0;
		ar >> m_F >> m_L >> m_eps >> m_atol;
	}
}
