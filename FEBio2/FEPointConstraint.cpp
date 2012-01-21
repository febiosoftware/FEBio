#include "stdafx.h"
#include "FEPointConstraint.h"
#include "fem.h"
#include "FESolidSolver.h"
#include "FEAnalysisStep.h"

//-----------------------------------------------------------------------------
FEPointConstraint::FEPointConstraint(FEM* pfem)
{
	m_pfem = pfem;
	m_node = -1;
	m_pel = 0;
}

//-----------------------------------------------------------------------------
void FEPointConstraint::Init()
{
	assert(m_node != -1);
	FEMesh& m = m_pfem->m_mesh;

	// get the nodal position in the reference state
	vec3d r = m.Node(m_node).m_r0;

	// find the element in which this node lies
	m_pel = m.FindSolidElement(r, m_rs);
	assert(m_pel);
}

//-----------------------------------------------------------------------------
void FEPointConstraint::Residual(vector<double> &R)
{
	int i;
	FEMesh& m = m_pfem->m_mesh;
	FEAnalysisStep* pstep = dynamic_cast<FEAnalysisStep*>(m_pfem->GetCurrentStep());
	FESolidSolver* psolver = dynamic_cast<FESolidSolver*>(pstep->m_psolver);

	// calculate H matrix
	double H[9], *r = m_rs;
	H[0] = 1.0;
	H[1] = -0.125*(1 - r[0])*(1 - r[1])*(1 - r[2]);
	H[2] = -0.125*(1 + r[0])*(1 - r[1])*(1 - r[2]);
	H[3] = -0.125*(1 + r[0])*(1 + r[1])*(1 - r[2]);
	H[4] = -0.125*(1 - r[0])*(1 + r[1])*(1 - r[2]);
	H[5] = -0.125*(1 - r[0])*(1 - r[1])*(1 + r[2]);
	H[6] = -0.125*(1 + r[0])*(1 - r[1])*(1 + r[2]);
	H[7] = -0.125*(1 + r[0])*(1 + r[1])*(1 + r[2]);
	H[8] = -0.125*(1 - r[0])*(1 + r[1])*(1 + r[2]);

	// get the nodal position
	vec3d x[9];
	x[0] = m.Node(m_node).m_rt;
	for (i=0; i<8; ++i) x[i+1] = m.Node(m_pel->m_node[i]).m_rt;

	// calculate the constraint
	vec3d c(0,0,0);
	for (i=0; i<9; ++i) c += x[i]*H[i];

	// calculate the force
	vec3d T = c*m_eps;

	// setup the LM matrix
	vector<int> LM(3*9), en(9);
	en[0] = m_node;
	LM[0] = m.Node(m_node).m_ID[DOF_X];
	LM[1] = m.Node(m_node).m_ID[DOF_Y];
	LM[2] = m.Node(m_node).m_ID[DOF_Z];
	for (i=0; i<8; ++i)
	{
		en[i+1] = m_pel->m_node[i];
		FENode& node = m.Node(en[i+1]);
		LM[(i+1)*3  ] = node.m_ID[DOF_X];
		LM[(i+1)*3+1] = node.m_ID[DOF_Y];
		LM[(i+1)*3+2] = node.m_ID[DOF_Z];
	}

	// set up nodal force vector
	vector<double> fe(3*9);
	for (int i=0; i<9; ++i)
	{
		fe[3*i  ] = -T.x*H[i];
		fe[3*i+1] = -T.y*H[i];
		fe[3*i+2] = -T.z*H[i];
	}

	// assemble residual
	psolver->AssembleResidual(en, LM, fe, R);
}

//-----------------------------------------------------------------------------
void FEPointConstraint::Stiffness()
{
	int i, j;
	FEMesh& m = m_pfem->m_mesh;
	FEAnalysisStep* pstep = dynamic_cast<FEAnalysisStep*>(m_pfem->GetCurrentStep());
	FESolidSolver* psolver = dynamic_cast<FESolidSolver*>(pstep->m_psolver);

	// calculate H matrix
	double H[9], *r = m_rs;
	H[0] = 1.0;
	H[1] = -0.125*(1 - r[0])*(1 - r[1])*(1 - r[2]);
	H[2] = -0.125*(1 + r[0])*(1 - r[1])*(1 - r[2]);
	H[3] = -0.125*(1 + r[0])*(1 + r[1])*(1 - r[2]);
	H[4] = -0.125*(1 - r[0])*(1 + r[1])*(1 - r[2]);
	H[5] = -0.125*(1 - r[0])*(1 - r[1])*(1 + r[2]);
	H[6] = -0.125*(1 + r[0])*(1 - r[1])*(1 + r[2]);
	H[7] = -0.125*(1 + r[0])*(1 + r[1])*(1 + r[2]);
	H[8] = -0.125*(1 - r[0])*(1 + r[1])*(1 + r[2]);


	// setup the LM matrix
	vector<int> LM(3*9), en(9);
	en[0] = m_node;
	LM[0] = m.Node(m_node).m_ID[DOF_X];
	LM[1] = m.Node(m_node).m_ID[DOF_Y];
	LM[2] = m.Node(m_node).m_ID[DOF_Z];
	for (i=0; i<8; ++i)
	{
		en[i+1] = m_pel->m_node[i];
		FENode& node = m.Node(en[i+1]);
		LM[(i+1)*3  ] = node.m_ID[DOF_X];
		LM[(i+1)*3+1] = node.m_ID[DOF_Y];
		LM[(i+1)*3+2] = node.m_ID[DOF_Z];
	}

	// setup stiffness matrix
	int ndof = 3*9;
	matrix ke(ndof, ndof); ke.zero();
	for (i=0; i<9; ++i)
		for (j=0; j<9; ++j)
		{
			ke[3*i  ][3*j  ] = m_eps*H[i]*H[j];
			ke[3*i+1][3*j+1] = m_eps*H[i]*H[j];
			ke[3*i+2][3*j+2] = m_eps*H[i]*H[j];
		}

	// assemble stiffness matrix
	psolver->AssembleStiffness(en, LM, ke);
}
