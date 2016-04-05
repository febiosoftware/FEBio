#include "stdafx.h"
#include "FEMicroMaterial2O.h"
#include "FECore/log.h"
#include "FESolidSolver2.h"
#include "FEElasticSolidDomain.h"
#include "FECore/FEAnalysis.h"
#include "FEBioXML/FEBioImport.h"
#include "FEBioPlot/FEBioPlotFile.h"
#include "FECore/tens3d.h"

//-----------------------------------------------------------------------------
FEMicroMaterialPoint2O::FEMicroMaterialPoint2O(FEMaterialPoint* mp) : FEMaterialPoint(mp)
{
	m_tau.zero();
	m_G.zero();

	m_inf_str_grad.zero();
	m_E.zero();
	m_H.zero();
	m_e.zero();
	m_h.zero();
	m_macro_energy = 0.;
	m_micro_energy = 0.;
	m_energy_diff = 0.;
	
	m_macro_energy_inc = 0.;
	m_micro_energy_inc = 0.;
	
	m_Ca.zero();
	m_Da.zero();
	m_Ea.zero();
	
	m_G_prev.zero();
}

//-----------------------------------------------------------------------------
//! Initialize material point data
void FEMicroMaterialPoint2O::Init(bool bflag)
{
	FEMaterialPoint::Init(bflag);
}

//-----------------------------------------------------------------------------
//! create a shallow copy
FEMaterialPoint* FEMicroMaterialPoint2O::Copy()
{
	FEMicroMaterialPoint2O* mmpt2O = new FEMicroMaterialPoint2O(m_pNext?m_pNext->Copy():0);
	mmpt2O->m_tau = m_tau;
	mmpt2O->m_G = m_G;
	mmpt2O->m_Ca = m_Ca;
	mmpt2O->m_Da = m_Da;
	mmpt2O->m_Ea = m_Ea;
	return mmpt2O;
}

//-----------------------------------------------------------------------------
//! serialize material point data
void FEMicroMaterialPoint2O::Serialize(DumpStream& ar)
{
	FEMaterialPoint::Serialize(ar);
	if (ar.IsSaving())
	{
		ar << m_Ca;
		ar << m_Da;
		ar << m_Ea;
	}
	else
	{
		ar >> m_Ca;
		ar >> m_Da;
		ar >> m_Ea;
	}
}

//=============================================================================

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_PARAMETER_LIST(FEMicroMaterial2O, FEElasticMaterial)
	ADD_PARAMETER(m_szrve    , FE_PARAM_STRING, "RVE"     );
	ADD_PARAMETER(m_szbc     , FE_PARAM_STRING, "bc_set"  );
	ADD_PARAMETER(m_bperiodic, FE_PARAM_BOOL  , "periodic");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEMicroMaterial2O::FEMicroMaterial2O(FEModel* pfem) : FEElasticMaterial2O(pfem)
{
	// initialize parameters
	m_szrve[0] = 0;
	m_szbc[0] = 0;
	m_bperiodic = false;

	AddProperty(&m_probe, "probe", false);
}

//-----------------------------------------------------------------------------
FEMicroMaterial2O::~FEMicroMaterial2O(void)
{
}

//-----------------------------------------------------------------------------
FEMaterialPoint* FEMicroMaterial2O::CreateMaterialPointData()
{
	return new FEMicroMaterialPoint2O(new FEElasticMaterialPoint);
}

//-----------------------------------------------------------------------------
bool FEMicroMaterial2O::Init()
{
	// initialize base class first
	if (FEElasticMaterial::Init() == false) return false;

	// load the master RVE model
	FEBioImport fim;
	if (fim.Load(m_mrve, m_szrve) == false)
	{
		return MaterialError("An error occured trying to read the RVE model from file %s.", m_szrve);
	}

	// the logfile is a shared resource between the master FEM and the RVE
	// in order not to corrupt the logfile we don't print anything for
	// the RVE problem.
	Logfile::MODE nmode = felog.GetMode();
	felog.SetMode(Logfile::NEVER);

	// initialize RVE
	if (m_mrve.InitRVE(m_bperiodic, m_szbc) == false) return MaterialError("An error occurred preparing RVE model");

	// reset the logfile mode
	felog.SetMode(nmode);

	return true;
}

//-----------------------------------------------------------------------------
//! Assign the prescribed displacement to the boundary nodes.
void FEMicroMaterial2O::UpdateBC(FEModel& rve, mat3d& F, tens3drs& G)
{
	// get the mesh
	FEMesh& m = rve.GetMesh();

	// assign new DC's for the boundary nodes
	FEPrescribedBC& dx = *rve.PrescribedBC(0);
	FEPrescribedBC& dy = *rve.PrescribedBC(1);
	FEPrescribedBC& dz = *rve.PrescribedBC(2);

	for (int i=0; i<(int) dx.Items(); ++i)
	{
		FENode& node = m.Node(dx.NodeID(i));
		const vec3d& r0 = node.m_r0;
		
		// Apply the second order boundary conditions to the RVE problem
		vec3d r1 = F*r0 + G.contractdyad1(r0)*0.5;

		// set the node scale
		dx.SetNodeScale(i, r1.x - r0.x);
		dy.SetNodeScale(i, r1.y - r0.y);
		dz.SetNodeScale(i, r1.z - r0.z);
	}


	if (m_bperiodic)
	{
		// get the "displacement" component of the deformation gradient
		mat3d U = F - mat3dd(1);

		// set the offset for the periodic BC's
		vec3d r[FEElement::MAX_NODES];

		// loop over periodic boundaries
		for (int i=0; i<3; ++i)
		{
			FEPeriodicBoundary2O* pc = dynamic_cast<FEPeriodicBoundary2O*>(rve.SurfacePairInteraction(i));
			assert(pc);
			pc->m_Fmacro = F;
			pc->m_Gmacro = G;
		}
	}
}

//-----------------------------------------------------------------------------
// The stiffness is evaluated at the same time the stress is evaluated so we 
// can just return it here. Note that this assumes that the stress function 
// is always called prior to the tangent function.
// LTE - Note that this function is not used in the second-order implemenetation
tens4ds FEMicroMaterial2O::Tangent(FEMaterialPoint &mp)
{
	FEMicroMaterialPoint2O& mmpt2O = *mp.ExtractData<FEMicroMaterialPoint2O>();
	return mmpt2O.m_Ca;
}

//-----------------------------------------------------------------------------
// LTE - Note that this function is not used in the second-order implemenetation
mat3ds FEMicroMaterial2O::Stress(FEMaterialPoint &mp)
{
	mat3ds sa; sa.zero();
	
	return sa;
}

//-----------------------------------------------------------------------------
//! Calculate the average stress from the RVE solution.
// LTE - Note that this function is not used in the second-order implemenetation
mat3ds FEMicroMaterial2O::AveragedStress(FEModel& rve, FEMaterialPoint &mp)
{
	mat3ds sa; sa.zero();
	
	return sa;
}

//-----------------------------------------------------------------------------
// The stiffness is evaluated at the same time the stress is evaluated so we 
// can just return it here. Note that this assumes that the stress function 
// is always called prior to the tangent function.
void  FEMicroMaterial2O::Tangent2O(FEMaterialPoint &mp, tens4ds& c, tens5ds& d, tens6ds& e)
{
	FEMicroMaterialPoint2O& mmpt = *mp.ExtractData<FEMicroMaterialPoint2O>();
	c = mmpt.m_Ca;
	d = mmpt.m_Da;
	e = mmpt.m_Ea;
}

//-----------------------------------------------------------------------------
void FEMicroMaterial2O::Stress2O(FEMaterialPoint &mp)
{
	// get the deformation gradient
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	FEMicroMaterialPoint2O& mmpt2O = *mp.ExtractData<FEMicroMaterialPoint2O>();
	mat3d F = pt.m_F;
	tens3drs G = mmpt2O.m_G;

	// apply the BC's
	UpdateBC(mmpt2O.m_rve, F, G);

	// solve the RVE
	bool bret = mmpt2O.m_rve.Solve();

	// make sure it converged
	if (bret == false) throw FEMultiScaleException();

	// calculate the averaged Cauchy stress
	AveragedStress2O(mmpt2O.m_rve, mp, pt.m_s, mmpt2O.m_tau);
	
	// calculate the averaged PK1 stress
//	AveragedStress2OPK1(mmpt2O.m_rve, mp, mmpt2O.m_PK1, mmpt2O.m_QK1);
	
	// calculate the averaged PK2 stress
	AveragedStress2OPK2(mmpt2O.m_rve, mp, mmpt2O.m_S, mmpt2O.m_T);

	// calculate the averaged stiffness
	AveragedStiffness(mmpt2O.m_rve, mp, mmpt2O.m_Ca, mmpt2O.m_Da, mmpt2O.m_Ea);

	// calculate the difference between the macro and micro energy for Hill-Mandel condition
	calc_energy_diff(mmpt2O.m_rve, mp);	
}

//-----------------------------------------------------------------------------
//! Calculate the average stress from the RVE solution.
void FEMicroMaterial2O::AveragedStress2O(FEModel& rve, FEMaterialPoint &mp, mat3ds &sa, tens3ds &taua)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	mat3d F = pt.m_F;
	double J = pt.m_J;

	// get the RVE mesh
	FEMesh& m = rve.GetMesh();

	mat3d s; s.zero();
	tens3ds tau; tau.zero();

	// for periodic BC's we take the reaction forces directly from the periodic constraints
	if (m_bperiodic)
	{
		// get the reaction for from the periodic constraints
		for (int i=0; i<3; ++i)
		{
			FEPeriodicBoundary2O* pbc = dynamic_cast<FEPeriodicBoundary2O*>(rve.SurfacePairInteraction(i));
			assert(pbc);
			FEPeriodicSurface& ss = pbc->m_ss;
			int N = ss.Nodes();
			for (int i=0; i<N; ++i)
			{
				FENode& node = ss.Node(i);
				vec3d f = ss.m_Fr[i];
				
				// We multiply by two since the reaction forces are only stored at the slave surface 
				// and we also need to sum over the master nodes (NOTE: should I figure out a way to 
				// store the reaction forces on the master nodes as well?)
				s += (f & node.m_rt)*2.0;

				vec3d x; x = node.m_rt;
		
				tau += dyad3s(x, f, x)*2.0;
			}
		}
	}

	// get the reaction force vector from the solid solver
	// (We also need to do this for the periodic BC, since at the prescribed nodes,
	// the contact forces will be zero).
	FEPrescribedBC& dx = *rve.PrescribedBC(0);
	FEPrescribedBC& dy = *rve.PrescribedBC(1);
	FEPrescribedBC& dz = *rve.PrescribedBC(2);

	const int dof_X = rve.GetDOFIndex("x");
	const int dof_Y = rve.GetDOFIndex("y");
	const int dof_Z = rve.GetDOFIndex("z");
	FEAnalysis* pstep = rve.GetCurrentStep();
	FESolidSolver2* ps = dynamic_cast<FESolidSolver2*>(pstep->GetFESolver());
	assert(ps);
	vector<double>& R = ps->m_Fr;

	int N = dx.Items();
	for (int i=0; i<N; ++i)
	{
		FENode& n = m.Node(dx.NodeID(i));
		vec3d f;
		f.x = R[-n.m_ID[dof_X]-2];
		f.y = R[-n.m_ID[dof_Y]-2];
		f.z = R[-n.m_ID[dof_Z]-2];
				
		s += (f & n.m_rt);

		vec3d x = n.m_rt;
		
		tau += dyad3s(x, f, x);
	}

	double V0 = m_mrve.InitialVolume();
	sa = s.sym() / (J*V0);
	taua = tau / (2*J*V0);
}

//-----------------------------------------------------------------------------
//! Calculate the average stiffness from the RVE solution.
void FEMicroMaterial2O::AveragedStiffness(FEModel& rve, FEMaterialPoint &mp, tens4ds& c, tens5ds& d, tens6ds& e)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// get the mesh
	FEMesh& m = rve.GetMesh();

	// get the solver
	FEAnalysis* pstep = rve.GetCurrentStep();
	FESolidSolver2* ps = dynamic_cast<FESolidSolver2*>(pstep->GetFESolver());

	// the element's stiffness matrix
	matrix ke;

	// element's residual
	vector<double> fe;

	// get deformation gradient and its inverse
	mat3d F = pt.m_F;
	mat3d Fi = F.inverse();
	double J = pt.m_J;

	// get the stress
	//mat3ds s = pt.m_s;

	// calculate the center point
	vec3d rc(0,0,0);
	for (int k=0; k<m.Nodes(); ++k) rc += m.Node(k).m_rt;
	rc /= (double) m.Nodes();

	// LTE - Calculate the initial center point
	vec3d rc0(0,0,0);
	for (int k=0; k<m.Nodes(); ++k) rc0 += m.Node(k).m_r0;
	rc0 /= (double) m.Nodes();

	c.zero();
	d.zero();
	e.zero();

	// LTE - elasticity tensor
	double D[6][6] = {0};
		
	// calculate the stiffness matrix and residual
	for (int k=0; k<m.Domains(); ++k)
	{
		FEElasticSolidDomain& bd = dynamic_cast<FEElasticSolidDomain&>(m.Domain(k));
		int NS = bd.Elements();
		for (int n=0; n<NS; ++n)
		{
			FESolidElement& el = bd.Element(n);

			// create the element's stiffness matrix
			int ne = el.Nodes();
			int ndof = 3*ne;
			ke.resize(ndof, ndof);
			ke.zero();

			// calculate the element's stiffness matrix
			bd.ElementStiffness(rve, n, ke);

			// create the element's residual
			fe.assign(ndof, 0);

			// calculate the element's residual
			bd.ElementInternalForce(el, fe);

			// loop over the element's nodes
			for (int i=0; i<ne; ++i)
			{
				FENode& ni = m.Node(el.m_node[i]);
				for (int j=0; j<ne; ++j)
				{
					FENode& nj = m.Node(el.m_node[j]);
					if (m_mrve.IsBoundaryNode(el.m_node[i]) && m_mrve.IsBoundaryNode(el.m_node[j]))
					{
						// both nodes are boundary nodes
						// so grab the element's submatrix
						double K[3][3];
						K[0][0] = ke[3*i  ][3*j  ]; K[0][1] = ke[3*i  ][3*j+1]; K[0][2] = ke[3*i  ][3*j+2];
						K[1][0] = ke[3*i+1][3*j  ]; K[1][1] = ke[3*i+1][3*j+1]; K[1][2] = ke[3*i+1][3*j+2];
						K[2][0] = ke[3*i+2][3*j  ]; K[2][1] = ke[3*i+2][3*j+1]; K[2][2] = ke[3*i+2][3*j+2];

						// get the nodal positions
						vec3d ri = ni.m_rt;
						vec3d rj = nj.m_rt;
						
						double Ri[3] = { ri.x, ri.y, ri.z };
						double Rj[3] = { rj.x, rj.y, rj.z };

						// create the elasticity tensor
						D[0][0] += Ri[0]*K[0][0]*Rj[0]; 
						D[1][1] += Ri[1]*K[1][1]*Rj[1]; 
						D[2][2] += Ri[2]*K[2][2]*Rj[2]; 

						D[0][1] += Ri[0]*K[0][1]*Rj[1];
						D[0][2] += Ri[0]*K[0][2]*Rj[2];
						D[1][2] += Ri[1]*K[1][2]*Rj[2];

						D[0][3] += 0.5*(Ri[0]*K[0][0]*Rj[1] + Ri[0]*K[0][1]*Rj[0]);
						D[0][4] += 0.5*(Ri[0]*K[0][1]*Rj[2] + Ri[0]*K[0][2]*Rj[1]);
						D[0][5] += 0.5*(Ri[0]*K[0][0]*Rj[2] + Ri[0]*K[0][2]*Rj[0]);

						D[1][3] += 0.5*(Ri[1]*K[1][0]*Rj[1] + Ri[1]*K[1][1]*Rj[0]);
						D[1][4] += 0.5*(Ri[1]*K[1][1]*Rj[2] + Ri[1]*K[1][2]*Rj[1]);
						D[1][5] += 0.5*(Ri[1]*K[1][0]*Rj[2] + Ri[1]*K[1][2]*Rj[0]);

						D[2][3] += 0.5*(Ri[2]*K[2][0]*Rj[1] + Ri[2]*K[2][1]*Rj[0]);
						D[2][4] += 0.5*(Ri[2]*K[2][1]*Rj[2] + Ri[2]*K[2][2]*Rj[1]);
						D[2][5] += 0.5*(Ri[2]*K[2][0]*Rj[2] + Ri[2]*K[2][2]*Rj[0]);

						D[3][3] += 0.25*(Ri[0]*K[1][0]*Rj[1] + Ri[1]*K[0][0]*Rj[1] + Ri[0]*K[1][1]*Rj[0] + Ri[1]*K[0][1]*Rj[0]);
						D[3][4] += 0.25*(Ri[0]*K[1][1]*Rj[2] + Ri[1]*K[0][1]*Rj[2] + Ri[0]*K[1][2]*Rj[1] + Ri[1]*K[0][2]*Rj[1]);
						D[3][5] += 0.25*(Ri[0]*K[1][0]*Rj[2] + Ri[1]*K[0][0]*Rj[2] + Ri[0]*K[1][2]*Rj[0] + Ri[1]*K[0][2]*Rj[0]);

						D[4][4] += 0.25*(Ri[1]*K[2][1]*Rj[2] + Ri[2]*K[1][1]*Rj[2] + Ri[1]*K[2][2]*Rj[1] + Ri[2]*K[1][2]*Rj[1]);
						D[4][5] += 0.25*(Ri[1]*K[2][0]*Rj[2] + Ri[2]*K[1][0]*Rj[2] + Ri[1]*K[2][2]*Rj[0] + Ri[2]*K[1][2]*Rj[0]);

						D[5][5] += 0.25*(Ri[0]*K[2][0]*Rj[2] + Ri[2]*K[0][0]*Rj[2] + Ri[0]*K[2][2]*Rj[0] + Ri[2]*K[0][2]*Rj[0]);

						calculate_d2O(d, K, Ri, Rj);
						calculate_e2O(e, K, Ri, Rj);
					}
				}
				
			}
		}
	}

	// divide by volume
	double V0 = m_mrve.InitialVolume();
	c = tens4ds(D)/(pt.m_J * V0);
	d = d/(2.*pt.m_J * V0);
	e = e/(4.*pt.m_J * V0);
}

//-----------------------------------------------------------------------------
void FEMicroMaterial2O::calculate_d2O(tens5ds& d, double K[3][3], double Ri[3], double Rj[3] )
{
	d.d[ 0] += calc_5ds_comp(K, Ri, Rj, 1, 1, 1, 1, 1);
	d.d[ 1] += calc_5ds_comp(K, Ri, Rj, 1, 1, 1, 1, 2);
	d.d[ 2] += calc_5ds_comp(K, Ri, Rj, 1, 1, 1, 1, 3);
	d.d[ 3] += calc_5ds_comp(K, Ri, Rj, 1, 1, 1, 2, 2);
	d.d[ 4] += calc_5ds_comp(K, Ri, Rj, 1, 1, 1, 2, 3);
	d.d[ 5] += calc_5ds_comp(K, Ri, Rj, 1, 1, 1, 3, 3);
	d.d[ 6] += calc_5ds_comp(K, Ri, Rj, 1, 1, 2, 2, 2);
	d.d[ 7] += calc_5ds_comp(K, Ri, Rj, 1, 1, 2, 2, 3);
	d.d[ 8] += calc_5ds_comp(K, Ri, Rj, 1, 1, 2, 3, 3);
	d.d[ 9] += calc_5ds_comp(K, Ri, Rj, 1, 1, 3, 3, 3);
	
	d.d[10] += calc_5ds_comp(K, Ri, Rj, 1, 2, 1, 2, 2);
	d.d[11] += calc_5ds_comp(K, Ri, Rj, 1, 2, 1, 2, 3);
	d.d[12] += calc_5ds_comp(K, Ri, Rj, 1, 2, 1, 3, 3);
	d.d[13] += calc_5ds_comp(K, Ri, Rj, 1, 2, 2, 2, 2);
	d.d[14] += calc_5ds_comp(K, Ri, Rj, 1, 2, 2, 2, 3);
	d.d[15] += calc_5ds_comp(K, Ri, Rj, 1, 2, 2, 3, 3);
	d.d[16] += calc_5ds_comp(K, Ri, Rj, 1, 2, 3, 3, 3);

	d.d[17] += calc_5ds_comp(K, Ri, Rj, 1, 3, 1, 3, 3);
	d.d[18] += calc_5ds_comp(K, Ri, Rj, 1, 3, 2, 2, 2);
	d.d[19] += calc_5ds_comp(K, Ri, Rj, 1, 3, 2, 2, 3);
	d.d[20] += calc_5ds_comp(K, Ri, Rj, 1, 3, 2, 3, 3);
	d.d[21] += calc_5ds_comp(K, Ri, Rj, 1, 3, 3, 3, 3);

	d.d[22] += calc_5ds_comp(K, Ri, Rj, 2, 2, 2, 2, 2);
	d.d[23] += calc_5ds_comp(K, Ri, Rj, 2, 2, 2, 2, 3);
	d.d[24] += calc_5ds_comp(K, Ri, Rj, 2, 2, 2, 3, 3);
	d.d[25] += calc_5ds_comp(K, Ri, Rj, 2, 2, 3, 3, 3);

	d.d[26] += calc_5ds_comp(K, Ri, Rj, 2, 3, 2, 3, 3);

	d.d[27] += calc_5ds_comp(K, Ri, Rj, 2, 3, 3, 3, 3);

	d.d[28] += calc_5ds_comp(K, Ri, Rj, 3, 3, 3, 3, 3);
}

//-----------------------------------------------------------------------------
double FEMicroMaterial2O::calc_5ds_comp(double K[3][3], double Ri[3], double Rj[3], int i, int j, int k, int l, int m)
{
	i -= 1; j -= 1;  k -= 1; l -= 1; m -= 1; 
		
	return (1/2)*((1/24)*(   Ri[i]*K[j][k]*Rj[l]*Rj[m] + Ri[i]*K[j][m]*Rj[l]*Rj[k] + Ri[i]*K[j][l]*Rj[k]*Rj[m] + Ri[i]*K[j][k]*Rj[m]*Rj[l] + Ri[i]*K[j][l]*Rj[m]*Rj[k] + Ri[i]*K[j][m]*Rj[k]*Rj[l]
					       + Ri[j]*K[i][k]*Rj[l]*Rj[m] + Ri[j]*K[i][m]*Rj[l]*Rj[k] + Ri[j]*K[i][l]*Rj[k]*Rj[m] + Ri[j]*K[i][k]*Rj[m]*Rj[l] + Ri[j]*K[i][l]*Rj[m]*Rj[k] + Ri[j]*K[i][m]*Rj[k]*Rj[l]
					       + Ri[l]*K[m][i]*Rj[j]*Rj[k] + Ri[l]*K[m][k]*Rj[j]*Rj[i] + Ri[l]*K[m][j]*Rj[i]*Rj[k] + Ri[l]*K[m][i]*Rj[k]*Rj[j] + Ri[l]*K[m][j]*Rj[k]*Rj[i] + Ri[l]*K[m][k]*Rj[i]*Rj[j]
					       + Ri[m]*K[l][i]*Rj[j]*Rj[k] + Ri[m]*K[l][k]*Rj[j]*Rj[i] + Ri[m]*K[l][j]*Rj[i]*Rj[k] + Ri[m]*K[l][i]*Rj[k]*Rj[j] + Ri[m]*K[l][j]*Rj[k]*Rj[i] + Ri[m]*K[l][k]*Rj[i]*Rj[j])
		        + (1/24)*(   Ri[i]*Ri[j]*K[k][l]*Rj[m] + Ri[i]*Ri[j]*K[m][l]*Rj[k] + Ri[i]*Ri[j]*K[l][k]*Rj[m] + Ri[i]*Ri[j]*K[k][m]*Rj[l] + Ri[i]*Ri[j]*K[l][m]*Rj[k] + Ri[i]*Ri[j]*K[m][k]*Rj[l]
					       + Ri[j]*Ri[i]*K[k][l]*Rj[m] + Ri[j]*Ri[i]*K[m][l]*Rj[k] + Ri[j]*Ri[i]*K[l][k]*Rj[m] + Ri[j]*Ri[i]*K[k][m]*Rj[l] + Ri[j]*Ri[i]*K[l][m]*Rj[k] + Ri[j]*Ri[i]*K[m][k]*Rj[l]
					       + Ri[l]*Ri[m]*K[i][j]*Rj[k] + Ri[l]*Ri[m]*K[k][j]*Rj[i] + Ri[l]*Ri[m]*K[j][i]*Rj[k] + Ri[l]*Ri[m]*K[i][k]*Rj[j] + Ri[l]*Ri[m]*K[j][k]*Rj[i] + Ri[l]*Ri[m]*K[k][i]*Rj[j]
					       + Ri[m]*Ri[l]*K[i][j]*Rj[k] + Ri[m]*Ri[l]*K[k][j]*Rj[i] + Ri[m]*Ri[l]*K[j][i]*Rj[k] + Ri[m]*Ri[l]*K[i][k]*Rj[j] + Ri[m]*Ri[l]*K[j][k]*Rj[i] + Ri[m]*Ri[l]*K[k][i]*Rj[j]));
}

//-----------------------------------------------------------------------------
void FEMicroMaterial2O::calculate_e2O(tens6ds& e, double K[3][3], double Ri[3], double Rj[3] )
{
	e.d[ 0] += calc_6ds_comp(K, Ri, Rj, 1, 1, 1, 1, 1, 1);
	e.d[ 1] += calc_6ds_comp(K, Ri, Rj, 1, 1, 1, 1, 1, 2);
	e.d[ 2] += calc_6ds_comp(K, Ri, Rj, 1, 1, 1, 1, 1, 3);
	e.d[ 3] += calc_6ds_comp(K, Ri, Rj, 1, 1, 1, 1, 2, 2);
	e.d[ 4] += calc_6ds_comp(K, Ri, Rj, 1, 1, 1, 1, 2, 3);
	e.d[ 5] += calc_6ds_comp(K, Ri, Rj, 1, 1, 1, 1, 3, 3);
	e.d[ 6] += calc_6ds_comp(K, Ri, Rj, 1, 1, 1, 2, 2, 2);
	e.d[ 7] += calc_6ds_comp(K, Ri, Rj, 1, 1, 1, 2, 2, 3);
	e.d[ 8] += calc_6ds_comp(K, Ri, Rj, 1, 1, 1, 2, 3, 3);
	e.d [9] += calc_6ds_comp(K, Ri, Rj, 1, 1, 1, 3, 3, 3);
	
	e.d[10] += calc_6ds_comp(K, Ri, Rj, 1, 1, 2, 1, 2, 2);
	e.d[11] += calc_6ds_comp(K, Ri, Rj, 1, 1, 2, 1, 2, 3);
	e.d[12] += calc_6ds_comp(K, Ri, Rj, 1, 1, 2, 1, 3, 3);
	e.d[13] += calc_6ds_comp(K, Ri, Rj, 1, 1, 2, 2, 2, 2);
	e.d[14] += calc_6ds_comp(K, Ri, Rj, 1, 1, 2, 2, 2, 3);
	e.d[15] += calc_6ds_comp(K, Ri, Rj, 1, 1, 2, 2, 3, 3);
	e.d[16] += calc_6ds_comp(K, Ri, Rj, 1, 1, 2, 3, 3, 3);

	e.d[17] += calc_6ds_comp(K, Ri, Rj, 1, 1, 3, 1, 2, 2);
	e.d[18] += calc_6ds_comp(K, Ri, Rj, 1, 1, 3, 1, 2, 3);
	e.d[19] += calc_6ds_comp(K, Ri, Rj, 1, 1, 3, 1, 3, 3);
	e.d[20] += calc_6ds_comp(K, Ri, Rj, 1, 1, 3, 2, 2, 2);
	e.d[21] += calc_6ds_comp(K, Ri, Rj, 1, 1, 3, 2, 2, 3);
	e.d[22] += calc_6ds_comp(K, Ri, Rj, 1, 1, 3, 2, 3, 3);
	e.d[23] += calc_6ds_comp(K, Ri, Rj, 1, 1, 3, 3, 3, 3);

	e.d[24] += calc_6ds_comp(K, Ri, Rj, 1, 2, 2, 1, 3, 3);
	e.d[25] += calc_6ds_comp(K, Ri, Rj, 1, 2, 2, 2, 2, 2);
	e.d[26] += calc_6ds_comp(K, Ri, Rj, 1, 2, 2, 2, 2, 3);
	e.d[27] += calc_6ds_comp(K, Ri, Rj, 1, 2, 2, 2, 3, 3);
	e.d[28] += calc_6ds_comp(K, Ri, Rj, 1, 2, 2, 3, 3, 3);

	e.d[29] += calc_6ds_comp(K, Ri, Rj, 1, 2, 3, 1, 3, 3);
	e.d[30] += calc_6ds_comp(K, Ri, Rj, 1, 2, 3, 2, 2, 2);
	e.d[31] += calc_6ds_comp(K, Ri, Rj, 1, 2, 3, 2, 2, 3);
	e.d[32] += calc_6ds_comp(K, Ri, Rj, 1, 2, 3, 2, 3, 3);
	e.d[33] += calc_6ds_comp(K, Ri, Rj, 1, 2, 3, 3, 3, 3);

	e.d[34] += calc_6ds_comp(K, Ri, Rj, 1, 3, 3, 2, 2, 2);
	e.d[35] += calc_6ds_comp(K, Ri, Rj, 1, 3, 3, 2, 2, 3);
	e.d[36] += calc_6ds_comp(K, Ri, Rj, 1, 3, 3, 2, 3, 3);
	e.d[37] += calc_6ds_comp(K, Ri, Rj, 1, 3, 3, 3, 3, 3);

	e.d[38] += calc_6ds_comp(K, Ri, Rj, 2, 2, 2, 2, 2, 2);
	e.d[39] += calc_6ds_comp(K, Ri, Rj, 2, 2, 2, 2, 2, 3);
	e.d[40] += calc_6ds_comp(K, Ri, Rj, 2, 2, 2, 2, 3, 3);
	e.d[41] += calc_6ds_comp(K, Ri, Rj, 2, 2, 2, 3, 3, 3);

	e.d[42] += calc_6ds_comp(K, Ri, Rj, 2, 2, 3, 2, 3, 3);
	e.d[43] += calc_6ds_comp(K, Ri, Rj, 2, 2, 3, 3, 3, 3);

	e.d[44] += calc_6ds_comp(K, Ri, Rj, 2, 3, 3, 3, 3, 3);

	e.d[45] += calc_6ds_comp(K, Ri, Rj, 3, 3, 3, 3, 3, 3);
}

//-----------------------------------------------------------------------------
double FEMicroMaterial2O::calc_6ds_comp(double K[3][3], double Ri[3], double Rj[3], int i, int j, int k, int l, int m, int n)
{
	i -= 1; j -= 1;  k -= 1; l -= 1; m -= 1; n -= 1;
		
	return (1/72)*(   Ri[i]*Ri[j]*K[k][l]*Rj[m]*Rj[n] + Ri[i]*Ri[j]*K[k][n]*Rj[m]*Rj[l] + Ri[i]*Ri[j]*K[k][m]*Rj[l]*Rj[n] + Ri[i]*Ri[j]*K[k][l]*Rj[n]*Rj[m] + Ri[i]*Ri[j]*K[k][m]*Rj[n]*Rj[l] + Ri[i]*Ri[j]*K[k][n]*Rj[l]*Rj[m]     
	                + Ri[k]*Ri[j]*K[i][l]*Rj[m]*Rj[n] + Ri[k]*Ri[j]*K[i][n]*Rj[m]*Rj[l] + Ri[k]*Ri[j]*K[i][m]*Rj[l]*Rj[n] + Ri[k]*Ri[j]*K[i][l]*Rj[n]*Rj[m] + Ri[k]*Ri[j]*K[i][m]*Rj[n]*Rj[l] + Ri[k]*Ri[j]*K[i][n]*Rj[l]*Rj[m]
				    + Ri[j]*Ri[i]*K[k][l]*Rj[m]*Rj[n] + Ri[j]*Ri[i]*K[k][n]*Rj[m]*Rj[l] + Ri[j]*Ri[i]*K[k][m]*Rj[l]*Rj[n] + Ri[j]*Ri[i]*K[k][l]*Rj[n]*Rj[m] + Ri[j]*Ri[i]*K[k][m]*Rj[n]*Rj[l] + Ri[j]*Ri[i]*K[k][n]*Rj[l]*Rj[m]
					+ Ri[i]*Ri[k]*K[j][l]*Rj[m]*Rj[n] + Ri[i]*Ri[k]*K[j][n]*Rj[m]*Rj[l] + Ri[i]*Ri[k]*K[j][m]*Rj[l]*Rj[n] + Ri[i]*Ri[k]*K[j][l]*Rj[n]*Rj[m] + Ri[i]*Ri[k]*K[j][m]*Rj[n]*Rj[l] + Ri[i]*Ri[k]*K[j][n]*Rj[l]*Rj[m]
					+ Ri[j]*Ri[k]*K[i][l]*Rj[m]*Rj[n] + Ri[j]*Ri[k]*K[i][n]*Rj[m]*Rj[l] + Ri[j]*Ri[k]*K[i][m]*Rj[l]*Rj[n] + Ri[j]*Ri[k]*K[i][l]*Rj[n]*Rj[m] + Ri[j]*Ri[k]*K[i][m]*Rj[n]*Rj[l] + Ri[j]*Ri[k]*K[i][n]*Rj[l]*Rj[m]
					+ Ri[k]*Ri[i]*K[j][l]*Rj[m]*Rj[n] + Ri[k]*Ri[i]*K[j][n]*Rj[m]*Rj[l] + Ri[k]*Ri[i]*K[j][m]*Rj[l]*Rj[n] + Ri[k]*Ri[i]*K[j][l]*Rj[n]*Rj[m] + Ri[k]*Ri[i]*K[j][m]*Rj[n]*Rj[l] + Ri[k]*Ri[i]*K[j][n]*Rj[l]*Rj[m]
					+ Ri[l]*Ri[m]*K[n][i]*Rj[j]*Rj[k] + Ri[l]*Ri[m]*K[n][k]*Rj[j]*Rj[i] + Ri[l]*Ri[m]*K[n][j]*Rj[i]*Rj[k] + Ri[l]*Ri[m]*K[n][i]*Rj[k]*Rj[j] + Ri[l]*Ri[m]*K[n][j]*Rj[k]*Rj[i] + Ri[l]*Ri[m]*K[n][k]*Rj[i]*Rj[j]
					+ Ri[n]*Ri[m]*K[l][i]*Rj[j]*Rj[k] + Ri[n]*Ri[m]*K[l][k]*Rj[j]*Rj[i] + Ri[n]*Ri[m]*K[l][j]*Rj[i]*Rj[k] + Ri[n]*Ri[m]*K[l][i]*Rj[k]*Rj[j] + Ri[n]*Ri[m]*K[l][j]*Rj[k]*Rj[i] + Ri[n]*Ri[m]*K[l][k]*Rj[i]*Rj[j]
					+ Ri[m]*Ri[l]*K[n][i]*Rj[j]*Rj[k] + Ri[m]*Ri[l]*K[n][k]*Rj[j]*Rj[i] + Ri[m]*Ri[l]*K[n][j]*Rj[i]*Rj[k] + Ri[m]*Ri[l]*K[n][i]*Rj[k]*Rj[j] + Ri[m]*Ri[l]*K[n][j]*Rj[k]*Rj[i] + Ri[m]*Ri[l]*K[n][k]*Rj[i]*Rj[j]
					+ Ri[l]*Ri[n]*K[m][i]*Rj[j]*Rj[k] + Ri[l]*Ri[n]*K[m][k]*Rj[j]*Rj[i] + Ri[l]*Ri[n]*K[m][j]*Rj[i]*Rj[k] + Ri[l]*Ri[n]*K[m][i]*Rj[k]*Rj[j] + Ri[l]*Ri[n]*K[m][j]*Rj[k]*Rj[i] + Ri[l]*Ri[n]*K[m][k]*Rj[i]*Rj[j]
					+ Ri[m]*Ri[n]*K[l][i]*Rj[j]*Rj[k] + Ri[m]*Ri[n]*K[l][k]*Rj[j]*Rj[i] + Ri[m]*Ri[n]*K[l][j]*Rj[i]*Rj[k] + Ri[m]*Ri[n]*K[l][i]*Rj[k]*Rj[j] + Ri[m]*Ri[n]*K[l][j]*Rj[k]*Rj[i] + Ri[m]*Ri[n]*K[l][k]*Rj[i]*Rj[j]
					+ Ri[n]*Ri[l]*K[m][i]*Rj[j]*Rj[k] + Ri[n]*Ri[l]*K[m][k]*Rj[j]*Rj[i] + Ri[n]*Ri[l]*K[m][j]*Rj[i]*Rj[k] + Ri[n]*Ri[l]*K[m][i]*Rj[k]*Rj[j] + Ri[n]*Ri[l]*K[m][j]*Rj[k]*Rj[i] + Ri[n]*Ri[l]*K[m][k]*Rj[i]*Rj[j]);
}

//-----------------------------------------------------------------------------
//! Calculate the energy difference between the RVE problem and the macro material point
void FEMicroMaterial2O::calc_energy_diff(FEModel& rve, FEMaterialPoint& mp)
{
	// get the deformation gradient and deformation hessian
/*	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	FEMicroMaterialPoint2O& mmpt2O = *mp.ExtractData<FEMicroMaterialPoint2O>();
	mat3d F = pt.m_F;
	tens3drs G = mmpt2O.m_G;

	mat3d Ftrans = F.transpose();
	tens3dls Gtrans = G.transpose();

	mat3d Finv = F.inverse();
	mat3d Finvtrans = Finv.transpose();
	tens3drs Ginv = G; Ginv.contractleg2(Finv,1); Ginv.contractleg2(Finv,2); Ginv.contractleg2(Finv,3);
	tens3dls Ginvtrans = Ginv.transpose();
	
	// calculate infinitesimal strain
	tens3d inf_strain_grad_nosym;

	inf_strain_grad_nosym.d[0] =  G.d[0];
	inf_strain_grad_nosym.d[1] =  G.d[1];
	inf_strain_grad_nosym.d[2] =  G.d[2];
	inf_strain_grad_nosym.d[3] =  0.5*(G.d[1] + G.d[6]);
	inf_strain_grad_nosym.d[4] =  0.5*(G.d[3] + G.d[7]);
	inf_strain_grad_nosym.d[5] =  0.5*(G.d[4] + G.d[8]);
	inf_strain_grad_nosym.d[6] =  0.5*(G.d[2] + G.d[12]);
	inf_strain_grad_nosym.d[7] =  0.5*(G.d[4] + G.d[13]);
	inf_strain_grad_nosym.d[8] =  0.5*(G.d[5] + G.d[14]);
	
	inf_strain_grad_nosym.d[9] =  0.5*(G.d[6] + G.d[1]);
	inf_strain_grad_nosym.d[10] = 0.5*(G.d[7] + G.d[3]);
	inf_strain_grad_nosym.d[11] = 0.5*(G.d[8] + G.d[4]);
	inf_strain_grad_nosym.d[12] = G.d[7];
	inf_strain_grad_nosym.d[13] = G.d[9];
	inf_strain_grad_nosym.d[14] = G.d[10];
	inf_strain_grad_nosym.d[15] = 0.5*(G.d[8] + G.d[13]);
	inf_strain_grad_nosym.d[16] = 0.5*(G.d[10] + G.d[15]);
	inf_strain_grad_nosym.d[17] = 0.5*(G.d[11] + G.d[16]);
	
	inf_strain_grad_nosym.d[18] = 0.5*(G.d[12] + G.d[2]);
	inf_strain_grad_nosym.d[19] = 0.5*(G.d[13] + G.d[4]);
	inf_strain_grad_nosym.d[20] = 0.5*(G.d[14] + G.d[5]);
	inf_strain_grad_nosym.d[21] = 0.5*(G.d[13] + G.d[8]);
	inf_strain_grad_nosym.d[22] = 0.5*(G.d[15] + G.d[10]);
	inf_strain_grad_nosym.d[23] = 0.5*(G.d[16] + G.d[11]);
	inf_strain_grad_nosym.d[24] = G.d[14];
	inf_strain_grad_nosym.d[25] = G.d[16];
	inf_strain_grad_nosym.d[26] = G.d[17];

	mmpt2O.m_inf_str_grad = inf_strain_grad_nosym.symm();

	// calculate Green-Lagrange strain
	mmpt2O.m_E = ((Ftrans*F - mat3dd(1))*0.5).sym();
	mmpt2O.m_H = (((Gtrans*F) + (Ftrans*G))*0.5).symm();

	// calculate Euler-Almansi strain
	mmpt2O.m_e = ((mat3dd(1) - Finvtrans*Finv)*0.5).sym();
	mmpt2O.m_h = (((Ginvtrans*Finv) + (Finvtrans*Ginv))*-0.5).symm();
	
	// calculate the energy difference between macro point and RVE
	// to verify that we have satisfied the Hill-Mandel condition
	mat3d F_prev = mmpt2O.m_F_prev;
	tens3drs G_prev = mmpt2O.m_G_prev;

	// calculate the macroscopic strain energy increment according to PK1 stress
	mmpt2O.m_macro_energy_inc = mmpt2O.m_PK1.dotdot(pt.m_F - F_prev) + mmpt2O.m_QK1.tripledot(mmpt2O.m_G - mmpt2O.m_G_prev);

	//// calculate the macroscopic strain energy increment according to PK2 stress
	///*mat3ds E_prev = ((F_prev.transpose()*F_prev - mat3dd(1))*0.5).sym();
	//tens3ds H_prev = ((G_prev.transpose().multiply2right(F).LStoUnsym() + G_prev.multiply2left(Ftrans).RStoUnsym())*0.5).symm();
	//mmpt2O.m_macro_energy_inc = mmpt2O.m_S.dotdot(mmpt2O.m_E - E_prev) + mmpt2O.m_T.tripledot(mmpt2O.m_H - H_prev);
	//
	//// calculate the macroscopic strain energy increment according to Cauchy stress
	///*mat3d Finv_prev = F_prev.inverse();
	//mat3d Finvtrans_prev = Finv_prev.transpose();
	//mat3ds e_prev = ((mat3dd(1) - F_prev.transinv()*F_prev.inverse())*0.5).sym();
	//tens3drs Ginv_prev; Ginv_prev = G_prev; Ginv_prev.contractleg2(Finv_prev,1); Ginv_prev.contractleg2(Finv_prev,2); Ginv_prev.contractleg2(Finv_prev,3);
	//tens3dls Ginvtrans_prev = Ginv_prev.transpose();
	//tens3ds h_prev = ((Ginvtrans_prev.multiply2right(Finv_prev).LStoUnsym() + Ginv_prev.multiply2left(Finvtrans_prev).RStoUnsym())*-0.5).symm();
	//mmpt2O.m_macro_energy_inc = pt.m_s.dotdot(mmpt2O.m_e - e_prev) + mmpt2O.m_tau.tripledot(mmpt2O.m_h - h_prev);

	// calculate the microscopic strain energy increment
	double rve_energy_avg = 0.;
	int nint; 
	double* w;
	
	mat3d rve_F;		
	mat3d rve_F_prev;
	double J0 = 0.;
	double V0 = 0.;
	mat3d rve_PK1;
	
	mat3ds rve_S;
	mat3ds rve_E;
	mat3ds rve_E_prev;
		
	mat3ds rve_s;
	mat3ds rve_e;
	mat3ds rve_e_prev;
	double J = 0.;
	double v = 0.;

	FEMesh& m = mmpt2O.m_rve.GetMesh();

	for (int k=0; k<m.Domains(); ++k)
	{
		FESolidDomain& dom = static_cast<FESolidDomain&>(m.Domain(k));

		for (int i=0; i<dom.Elements(); ++i)
		{
			FESolidElement& el = dom.Element(i);
			
			nint = el.GaussPoints();
			w = el.GaussWeights();
			
			for (int n=0; n<nint; ++n)
			{
				FEMaterialPoint& mp = *el.GetMaterialPoint(n);
				FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
				FEMicroMaterialPoint2O& mmpt = *mp.ExtractData<FEMicroMaterialPoint2O>();

				rve_F = pt.m_F;
				rve_F_prev = mmpt.m_F_prev;
				rve_s = pt.m_s;

				//// calculate microscopic strain energy according to PK1 stress
				//rve_PK1 = rve_F.det()*rve_s*rve_F.transinv();
				//J0 = dom.detJ0(el, n);		
				//V0 += J0*w[n];
				//rve_energy_avg += rve_PK1.dotdot(rve_F - rve_F_prev)*J0*w[n];
				//				
				//// calculate microscopic strain energy according to PK2 stress
				///*rve_E = ((rve_F.transpose()*rve_F - mat3dd(1))*0.5).sym();
				//rve_E_prev = ((rve_F_prev.transpose()*rve_F_prev - mat3dd(1))*0.5).sym();
				//rve_PK1 = rve_F.det()*rve_s*rve_F.transinv();
				//rve_S = (rve_F.inverse()*rve_PK1).sym();
				//J0 = dom.detJ0(el, n);		
				//V0 += J0*w[n];
				//rve_energy_avg += rve_S.dotdot(rve_E - rve_E_prev)*J0*w[n];
				//
				//// calculate microscopic strain energy according to Cauchy stress
				///*rve_s = rve_pt.m_s;
				//rve_e = ((mat3dd(1) - rve_F.transinv()*rve_F.inverse())*0.5).sym();
				//rve_e_prev = ((mat3dd(1) - rve_F_prev.transinv()*rve_F_prev.inverse())*0.5).sym();
				//J = dom.detJt(el, n);		
				//v += J*w[n];		
				//rve_energy_avg += rve_s.dotdot(rve_e - rve_e_prev)*J*w[n];
			}
		}
	}

	mmpt2O.m_micro_energy_inc = rve_energy_avg/V0;
*/
}


//-----------------------------------------------------------------------------
//! Calculate the average stress from the RVE solution.
mat3d FEMicroMaterial2O::AveragedStressPK1(FEModel& rve, FEMaterialPoint &mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	mat3d F = pt.m_F;
	double J = pt.m_J;
	
	// get the RVE mesh
	FEMesh& m = rve.GetMesh();

	mat3d PK1; PK1.zero();

	// for periodic BC's we take the reaction forces directly from the periodic constraints
	if (m_bperiodic)
	{
		// get the reaction for from the periodic constraints
		for (int i=0; i<3; ++i)
		{
			FEPeriodicBoundary1O* pbc = dynamic_cast<FEPeriodicBoundary1O*>(rve.SurfacePairInteraction(i));
			assert(pbc);
			FEPeriodicSurface& ss = pbc->m_ss;
			int N = ss.Nodes();
			for (int i=0; i<N; ++i)
			{
				FENode& node = ss.Node(i);
				vec3d f = ss.m_Fr[i];

				// We multiply by two since the reaction forces are only stored at the slave surface 
				// and we also need to sum over the master nodes (NOTE: should I figure out a way to 
				// store the reaction forces on the master nodes as well?)
				PK1 += (f & node.m_r0)*2.0;
			}
		}
	}

	// get the reaction force vector from the solid solver
	// (We also need to do this for the periodic BC, since at the prescribed nodes,
	// the contact forces will be zero). 
	const int dof_X = rve.GetDOFIndex("x");
	const int dof_Y = rve.GetDOFIndex("y");
	const int dof_Z = rve.GetDOFIndex("z");
	FEAnalysis* pstep = rve.GetCurrentStep();
	FESolidSolver2* ps = dynamic_cast<FESolidSolver2*>(pstep->GetFESolver());
	assert(ps);
	vector<double>& R = ps->m_Fr;
	FEPrescribedBC& dc = *rve.PrescribedBC(0);
	int nitems = dc.Items();
	for (int i=0; i<nitems; ++i)
	{
		FENode& n = m.Node(dc.NodeID(i));
		vec3d f;
		f.x = R[-n.m_ID[dof_X]-2];
		f.y = R[-n.m_ID[dof_Y]-2];
		f.z = R[-n.m_ID[dof_Z]-2];
		PK1 += f & n.m_r0;
	}

	double V0 = m_mrve.InitialVolume();
	return PK1 / V0;
}

//-----------------------------------------------------------------------------
//! Calculate the average stress from the RVE solution.
void FEMicroMaterial2O::AveragedStress2OPK1(FEModel& rve, FEMaterialPoint &mp, mat3d &PK1a, tens3drs &QK1a)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	mat3d F = pt.m_F;
	double J = pt.m_J;

	// get the RVE mesh
	FEMesh& m = rve.GetMesh();

	mat3d PK1; PK1.zero();
	tens3drs QK1; QK1.zero();

	// for periodic BC's we take the reaction forces directly from the periodic constraints
	if (m_bperiodic)
	{
		// get the reaction for from the periodic constraints
		for (int i=0; i<3; ++i)
		{
			FEPeriodicBoundary2O* pbc = dynamic_cast<FEPeriodicBoundary2O*>(rve.SurfacePairInteraction(i));
			assert(pbc);
			FEPeriodicSurface& ss = pbc->m_ss;
			int N = ss.Nodes();
			for (int i=0; i<N; ++i)
			{
				FENode& node = ss.Node(i);
				vec3d f = ss.m_Fr[i];

				// We multiply by two since the reaction forces are only stored at the slave surface 
				// and we also need to sum over the master nodes (NOTE: should I figure out a way to 
				// store the reaction forces on the master nodes as well?)
				PK1 += (f & node.m_r0)*2.0;

				vec3d X = node.m_r0;
		
				QK1 += dyad3rs(f, X)*2.0;
			}
		}
	}

	// get the reaction force vector from the solid solver
	// (We also need to do this for the periodic BC, since at the prescribed nodes,
	// the contact forces will be zero). 
	FEPrescribedBC& dx = *rve.PrescribedBC(0);
	FEPrescribedBC& dy = *rve.PrescribedBC(1);
	FEPrescribedBC& dz = *rve.PrescribedBC(2);

	const int dof_X = rve.GetDOFIndex("x");
	const int dof_Y = rve.GetDOFIndex("y");
	const int dof_Z = rve.GetDOFIndex("z");
	FEAnalysis* pstep = rve.GetCurrentStep();
	FESolidSolver2* ps = dynamic_cast<FESolidSolver2*>(pstep->GetFESolver());
	assert(ps);
	vector<double>& R = ps->m_Fr;
	int N = dx.Items();
	for (int i=0; i<N; ++i)
	{
		FENode& n = m.Node(dx.NodeID(i));
		vec3d f;
		f.x = R[-n.m_ID[dof_X]-2];
		f.y = R[-n.m_ID[dof_Y]-2];
		f.z = R[-n.m_ID[dof_Z]-2];
		
		PK1 += f & n.m_r0;
		vec3d X; X = n.m_r0;
		
		QK1 += dyad3rs(f, X);
	}

	double V0 = m_mrve.InitialVolume();
	PK1a = PK1 / V0;
	QK1a = QK1 / (2*V0);
}

//-----------------------------------------------------------------------------
//! Calculate the average stress from the RVE solution.
void FEMicroMaterial2O::AveragedStress2OPK2(FEModel& rve, FEMaterialPoint &mp, mat3ds &Sa, tens3ds &Ta)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	mat3d F = pt.m_F;
	double J = pt.m_J;
	mat3d Finv = F.inverse();

	// get the RVE mesh
	FEMesh& m = rve.GetMesh();

	mat3d S; S.zero();
	tens3ds T; T.zero();

	// for periodic BC's we take the reaction forces directly from the periodic constraints
	if (m_bperiodic)
	{
		// get the reaction for from the periodic constraints
		for (int i=0; i<3; ++i)
		{
			FEPeriodicBoundary2O* pbc = dynamic_cast<FEPeriodicBoundary2O*>(rve.SurfacePairInteraction(i));
			assert(pbc);
			FEPeriodicSurface& ss = pbc->m_ss;
			int N = ss.Nodes();
			for (int i=0; i<N; ++i)
			{
				FENode& node = ss.Node(i);
				vec3d f = ss.m_Fr[i];
				vec3d f0 = Finv*f;

				// We multiply by two since the reaction forces are only stored at the slave surface 
				// and we also need to sum over the master nodes (NOTE: should I figure out a way to 
				// store the reaction forces on the master nodes as well?)
				S += (f0 & node.m_r0)*2.0;

				vec3d X = node.m_r0;
		
				T += dyad3s(X, f0, X)*2.0;
			}
		}
	}

	// get the reaction force vector from the solid solver
	// (We also need to do this for the periodic BC, since at the prescribed nodes,
	// the contact forces will be zero). 
	FEPrescribedBC& dx = *rve.PrescribedBC(0);
	FEPrescribedBC& dy = *rve.PrescribedBC(1);
	FEPrescribedBC& dz = *rve.PrescribedBC(2);

	const int dof_X = rve.GetDOFIndex("x");
	const int dof_Y = rve.GetDOFIndex("y");
	const int dof_Z = rve.GetDOFIndex("z");
	FEAnalysis* pstep = rve.GetCurrentStep();
	FESolidSolver2* ps = dynamic_cast<FESolidSolver2*>(pstep->GetFESolver());
	assert(ps);
	vector<double>& R = ps->m_Fr;

	int N = dx.Items();
	
	for (int i=0; i<N; ++i)
	{
		FENode& n = m.Node(dx.NodeID(i));
		vec3d f;
		f.x = R[-n.m_ID[dof_X]-2];
		f.y = R[-n.m_ID[dof_Y]-2];
		f.z = R[-n.m_ID[dof_Z]-2];
		vec3d f0 = Finv*f;
		
		S += f0 & n.m_r0;

		vec3d X = n.m_r0;
		
		T += dyad3s(X, f0, X);
	}

	double V0 = m_mrve.InitialVolume();
	Sa = S.sym() / V0;
	Ta = T / (2*V0);
}
