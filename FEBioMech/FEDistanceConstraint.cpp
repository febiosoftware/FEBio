#include "stdafx.h"
#include "FEDistanceConstraint.h"
#include <FECore/FEModel.h>
#include <FECore/log.h>

//-----------------------------------------------------------------------------
BEGIN_PARAMETER_LIST(FEDistanceConstraint, FENLConstraint);
	ADD_PARAMETER(m_blaugon, FE_PARAM_BOOL  , "laugon" ); 
	ADD_PARAMETER(m_atol   , FE_PARAM_DOUBLE, "augtol" );
	ADD_PARAMETER(m_eps    , FE_PARAM_DOUBLE, "penalty");
	ADD_PARAMETERV(m_node  , FE_PARAM_INTV  , 2, "node");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
//! constructor
FEDistanceConstraint::FEDistanceConstraint(FEModel* pfem) : FENLConstraint(pfem)
{
	m_eps = 0.0;
	m_atol = 0.01;
	m_blaugon = false;
	m_node[0] = -1;
	m_node[1] = -1;
	m_L0 = 0.0;
}

//-----------------------------------------------------------------------------
//! Initializes data structures. 
bool FEDistanceConstraint::Init()
{
	// get the FE mesh
	FEMesh& mesh = GetFEModel()->GetMesh();
	int NN = mesh.Nodes();

	// make sure the nodes are valid
	// (remember, they are one-based since they are defined in the input file)
	if ((m_node[0] <= 0)||(m_node[0] > NN)) return false;
	if ((m_node[1] <= 0)||(m_node[1] > NN)) return false;
	
	// get the initial position of the two nodes
	vec3d ra = mesh.Node(m_node[0] - 1).m_r0;
	vec3d rb = mesh.Node(m_node[1] - 1).m_r0;

	// set the initial length
	m_L0 = (ra - rb).norm();

	return true;
}

//-----------------------------------------------------------------------------
void FEDistanceConstraint::Residual(FEGlobalVector& R)
{
	// get the FE mesh
	FEMesh& mesh = GetFEModel()->GetMesh();

	// get the two nodes
	FENode& nodea = mesh.Node(m_node[0] - 1);
	FENode& nodeb = mesh.Node(m_node[1] - 1);

	// get the current position of the two nodes
	vec3d ra = nodea.m_rt;
	vec3d rb = nodeb.m_rt;

	// calculate the force
	double L = (ra - rb).norm();
	vec3d Fc = m_Fc + (ra - rb)*(m_eps*(L - m_L0));

	// setup the "element" force vector
	vector<double> fe(6);
	fe[0] = -Fc.x;
	fe[1] = -Fc.y;
	fe[2] = -Fc.z;
	fe[3] =  Fc.x;
	fe[4] =  Fc.y;
	fe[5] =  Fc.z;

	// setup the LM vector
	vector<int> lm(6);
	lm[0] = nodea.m_ID[DOF_X];
	lm[1] = nodea.m_ID[DOF_Y];
	lm[2] = nodea.m_ID[DOF_Z];
	lm[3] = nodeb.m_ID[DOF_X];
	lm[4] = nodeb.m_ID[DOF_Y];
	lm[5] = nodeb.m_ID[DOF_Z];

	// setup element vector
	vector<int> en(2);
	en[0] = m_node[0] - 1;
	en[1] = m_node[1] - 1;

	// add element force vector to global force vector
	R.Assemble(en, lm, fe);
}

//-----------------------------------------------------------------------------
void FEDistanceConstraint::StiffnessMatrix(FESolver* psolver)
{
	// get the FE mesh
	FEMesh& mesh = GetFEModel()->GetMesh();

	// get the two nodes
	FENode& nodea = mesh.Node(m_node[0] - 1);
	FENode& nodeb = mesh.Node(m_node[1] - 1);

	// get the current position of the two nodes
	vec3d ra = nodea.m_rt;
	vec3d rb = nodeb.m_rt;
	vec3d rab = ra - rb;

	// calculate the "gap"
	double L = rab.norm();
	double g = m_eps*(L - m_L0);

	// calculate the stiffness
	mat3d kab;
	kab[0][0] = (m_eps*rab.x*rab.x/L + g);
	kab[1][1] = (m_eps*rab.y*rab.y/L + g);
	kab[2][2] = (m_eps*rab.z*rab.z/L + g);
	kab[0][1] = kab[1][0] = (m_eps*rab.x*rab.y/L);
	kab[0][2] = kab[2][0] = (m_eps*rab.x*rab.z/L);
	kab[1][2] = kab[2][1] = (m_eps*rab.y*rab.z/L);

	// element stiffness matrix
	matrix ke(6, 6);
	ke.zero();
	ke[0][0] = kab[0][0]; ke[0][1] = kab[0][1]; ke[0][2] = kab[0][2]; ke[0][3] = -kab[0][0]; ke[0][4] = -kab[0][1]; ke[0][5] = -kab[0][2];
	ke[1][0] = kab[1][0]; ke[1][1] = kab[1][1]; ke[1][2] = kab[1][2]; ke[1][3] = -kab[1][0]; ke[1][4] = -kab[1][1]; ke[1][5] = -kab[1][2];
	ke[2][0] = kab[2][0]; ke[2][1] = kab[2][1]; ke[2][2] = kab[2][2]; ke[2][3] = -kab[2][0]; ke[2][4] = -kab[2][1]; ke[2][5] = -kab[2][2];

	ke[3][0] = -kab[0][0]; ke[3][1] = -kab[0][1]; ke[3][2] = -kab[0][2]; ke[3][3] = kab[0][0]; ke[3][4] = kab[0][1]; ke[3][5] = kab[0][2];
	ke[4][0] = -kab[1][0]; ke[4][1] = -kab[1][1]; ke[4][2] = -kab[1][2]; ke[4][3] = kab[1][0]; ke[4][4] = kab[1][1]; ke[4][5] = kab[1][2];
	ke[5][0] = -kab[2][0]; ke[5][1] = -kab[2][1]; ke[5][2] = -kab[2][2]; ke[5][3] = kab[2][0]; ke[5][4] = kab[2][1]; ke[5][5] = kab[2][2];

	// setup the LM vector
	vector<int> lm(6);
	lm[0] = nodea.m_ID[DOF_X];
	lm[1] = nodea.m_ID[DOF_Y];
	lm[2] = nodea.m_ID[DOF_Z];
	lm[3] = nodeb.m_ID[DOF_X];
	lm[4] = nodeb.m_ID[DOF_Y];
	lm[5] = nodeb.m_ID[DOF_Z];

	// setup element vector
	vector<int> en(2);
	en[0] = m_node[0] - 1;
	en[1] = m_node[1] - 1;

	// assemble element matrix in global stiffness matrix
	psolver->AssembleStiffness(en, lm, ke);
}

//-----------------------------------------------------------------------------
bool FEDistanceConstraint::Augment(int naug)
{
	// make sure we are augmenting
	if ((m_blaugon == false) || (m_atol <= 0.0)) return true;

	// get the FE mesh
	FEMesh& mesh = GetFEModel()->GetMesh();

	// get the two nodes
	FENode& nodea = mesh.Node(m_node[0] - 1);
	FENode& nodeb = mesh.Node(m_node[1] - 1);

	// get the current position of the two nodes
	vec3d ra = nodea.m_rt;
	vec3d rb = nodeb.m_rt;

	// calculate the force
	double L = (ra - rb).norm();
	vec3d Fc = m_Fc + (ra - rb)*(m_eps*(L - m_L0));

	double err = (m_Fc - Fc).norm() / Fc.norm();

	felog.printf("\ndistance constraint:\n");
	felog.printf("\tforce    = %lg, %lg, %lg\n", Fc.x, Fc.y, Fc.z);
	felog.printf("\terror    = %lg (%lg)\n", err, m_atol);
	felog.printf("\tdistance = %lg (%lg)\n", L, m_L0);

	// check convergence
	if (err < m_atol) return true;

	// update Lagrange multiplier
	m_Fc = Fc;

	return false;
}

//-----------------------------------------------------------------------------
void FEDistanceConstraint::Serialize(DumpFile& ar)
{
}

//-----------------------------------------------------------------------------
void FEDistanceConstraint::Reset()
{
	m_Fc = vec3d(0,0,0);
}

//-----------------------------------------------------------------------------
// This function is called when the FE model's state needs to be updated.
void FEDistanceConstraint::Update()
{

}
