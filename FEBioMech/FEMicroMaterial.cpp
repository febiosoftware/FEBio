#include "stdafx.h"
#include "FEMicroMaterial.h"
#include "FECore/FEElemElemList.h"
#include "FECore/log.h"
#include "FEBioMech/FESolidSolver.h"
#include "FEBioMech/FEElasticSolidDomain.h"
#include "FECore/FEAnalysis.h"
#include "FEBioXML/FEBioImport.h"

// define the material parameters
BEGIN_PARAMETER_LIST(FEMicroMaterial, FEElasticMaterial)
	ADD_PARAMETER(m_szrve, FE_PARAM_STRING, "RVE");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEMicroMaterial::FEMicroMaterial(FEModel* pfem) : FEElasticMaterial(pfem)
{
	m_szrve[0] = 0;
	m_brve = false;
}

//-----------------------------------------------------------------------------
FEMicroMaterial::~FEMicroMaterial(void)
{
}

//-----------------------------------------------------------------------------
void FEMicroMaterial::Init()
{
	// try to load the RVE model
	if (m_brve == false)
	{
		FEFEBioImport fim;
		if (fim.Load(m_rve, m_szrve) == false)
		{
			throw MaterialError("An error occured trying to read the RVE model from file %s.", m_szrve);
		}

		// make sure the RVE problem doesn't output anything to a plot file
		m_rve.GetCurrentStep()->SetPlotLevel(FE_PLOT_NEVER);

		// create the DC's for this RVE
		PrepRVE();

		m_brve = true;
	}
}

//-----------------------------------------------------------------------------
void FEMicroMaterial::PrepRVE()
{
	// first we need to find all the boundary nodes
	FEMesh& m = m_rve.GetMesh();
	int N = m.Nodes();
	vector<int> tag; tag.assign(N, 0);

	// create the element-element list
	FEElemElemList EEL;
	EEL.Create(&m);

	// use the E-E list to tag all exterior nodes
	int fn[FEElement::MAX_NODES], nf, M = 0;
	FEElasticSolidDomain& bd = dynamic_cast<FEElasticSolidDomain&>(m.Domain(0));
	for (int i=0; i<bd.Elements(); ++i, ++M)
	{
		FESolidElement& el = bd.Element(i);
		nf = m.Faces(el);
		for (int j=0; j<nf; ++j)
		{
			if (EEL.Neighbor(M, j) == 0)
			{
				// mark all nodes
				int nn = m.GetFace(el, j, fn);
				for (int k=0; k<nn; ++k) tag[fn[k]] = 1;
			}
		}
	}

	// count the nr of exterior nodes
	int NN = 0, i;
	for (i=0; i<N; ++i) if (tag[i] == 1) ++NN;

	assert(NN > 0);

	// create a load curve
	FELoadCurve* plc = new FELoadCurve;
	plc->SetInterpolation(FELoadCurve::LINEAR);
	plc->Add(0.0, 0.0);
	plc->Add(1.0, 1.0);
	m_rve.AddLoadCurve(plc);
	int NLC = m_rve.LoadCurves() - 1;

	// create the DC's
	NN = 0;
	m_rve.ClearBCs();
	for (i=0; i<N; ++i)
		if (tag[i] == 1)
		{
			for (int j=0; j<3; ++j, ++NN)
			{
				FEPrescribedBC* pdc = new FEPrescribedBC(&m_rve);
				pdc->bc = j;
				pdc->lc = NLC;
				pdc->node = i;
				pdc->s = 0;
				m_rve.AddPrescribedBC(pdc);
			}
		}

	// the logfile is a shared resource between the master FEM and the RVE
	// in order not to corrupt the logfile we don't print anything for
	// the RVE problem.
	Logfile::MODE nmode = felog.GetMode();
	felog.SetMode(Logfile::NEVER);

	// initialize RVE
	m_rve.Init();

	// calculate intial RVE volume
	m_V0 = 0;
	double ve;
	int nint, n;
	double* w, J;
	for (i=0; i<bd.Elements(); ++i)
	{
		FESolidElement& el = bd.Element(i);
		nint = el.GaussPoints();
		w = el.GaussWeights();
		ve = 0;
		for (n=0; n<nint; ++n)
		{
			FEElasticMaterialPoint& pt = *el.m_State[n]->ExtractData<FEElasticMaterialPoint>();
			J = bd.detJt(el, n);

			ve += J*w[n];
		}
		m_V0 += ve;
	}

	// reset the logfile mode
	felog.SetMode(nmode);
}

//-----------------------------------------------------------------------------
mat3ds FEMicroMaterial::Stress(FEMaterialPoint &mp)
{
	// get the deformation gradient
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	mat3d F = pt.m_F;

	// the logfile is a shared resource between the master FEM and the RVE
	// in order not to corrupt the logfile we don't print anything for
	// the RVE problem.
	Logfile::MODE nmode = felog.GetMode();
	felog.SetMode(Logfile::NEVER);

	// reset the RVE
	m_rve.Reset();

	// get the mesh
	FEMesh& m = m_rve.GetMesh();

	// assign new DC's for the boundary nodes
	int N = m_rve.PrescribedBCs()/3, i;
	for (i=0; i<N; ++i)
	{
		FEPrescribedBC& dx = *m_rve.PrescribedBC(3*i  );
		FEPrescribedBC& dy = *m_rve.PrescribedBC(3*i+1);
		FEPrescribedBC& dz = *m_rve.PrescribedBC(3*i+2);

		FENode& node = m.Node(dx.node);

		vec3d r0 = node.m_r0;
		vec3d r1 = F*r0;

		dx.s = r1.x - r0.x;
		dy.s = r1.y - r0.y;
		dz.s = r1.z - r0.z;
	}

	// solve the RVE
	bool bret = m_rve.Solve();

	// reset the logfile mode
	felog.SetMode(nmode);

	if (bret == false) throw FEMultiScaleException();

	// calculate the averaged stress
	return AveragedStress(pt);
}

//-----------------------------------------------------------------------------

mat3ds FEMicroMaterial::AveragedStress(FEMaterialPoint& mp)
{
	// get the deformation gradient
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	mat3d F = pt.m_F;
	double J = pt.m_J;

	// get the mesh
	FEMesh& m = m_rve.GetMesh();
/*
	mat3ds s(0);
	double V = 0, ve;
	int nint, n, i;
	double* w, J;
	FEElasticSolidDomain& bd = dynamic_cast<FEElasticSolidDomain&>(m.Domain(0));
	for (i=0; i<bd.Elements(); ++i)
	{
		FESolidElement& el = bd.Element(i);
		m.UnpackElement(el);
		nint = el.GaussPoints();
		w = el.GaussWeights();
		ve = 0;
		for (n=0; n<nint; ++n)
		{
			FEElasticMaterialPoint& pt = *el.m_State[n]->ExtractData<FEElasticMaterialPoint>();
			J = bd.detJt(el, n);

			ve += J*w[n];
			s += pt.s*(J*w[n]);
		}
		V += ve;
	}
	s /= V;
*/
	// get the reaction force vector from the solid solver
	FEAnalysis* pstep = m_rve.GetCurrentStep();
	FESolidSolver* ps = dynamic_cast<FESolidSolver*>(pstep->m_psolver);
	assert(ps);
	vector<double>& R = ps->m_Fr;
	mat3d T; T.zero();
	int nbc = m_rve.PrescribedBCs();
	for (int i=0; i<nbc/3; ++i)
	{
		FEPrescribedBC& dc = *m_rve.PrescribedBC(3*i);
		FENode& n = m.Node(dc.node);
		vec3d f;
		f.x = R[-n.m_ID[DOF_X]-2];
		f.y = R[-n.m_ID[DOF_Y]-2];
		f.z = R[-n.m_ID[DOF_Z]-2];
		T += f & n.m_rt;
	}
	mat3ds s = T.sym() / (J*m_V0);
	return s;
}

//-----------------------------------------------------------------------------
tens4ds FEMicroMaterial::Tangent(FEMaterialPoint &mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// get the mesh
	FEMesh& m = m_rve.GetMesh();

	// get the solver
	FEAnalysis* pstep = m_rve.GetCurrentStep();
	FESolidSolver* ps = dynamic_cast<FESolidSolver*>(pstep->m_psolver);

	// the element's stiffness matrix
	matrix ke;

	// element's residual
	vector<double> fe;

	// elasticity tensor
	double D[6][6] = {0};

	// get deformation gradient and its inverse
	mat3d F = pt.m_F;
	mat3d Fi = F.inverse();

	// get the stress
	mat3ds s = pt.m_s;

	// calculate the stiffness matrix and residual
	FEElasticSolidDomain& bd = dynamic_cast<FEElasticSolidDomain&>(m.Domain(0));
	int NS = bd.Elements(), i, j;
	for (int n=0; n<NS; ++n)
	{
		FESolidElement& e = bd.Element(n);

		// create the element's stiffness matrix
		int ne = e.Nodes();
		int ndof = 3*ne;
		ke.resize(ndof, ndof);
		ke.zero();

		// calculate the element's stiffness matrix
		bd.ElementStiffness(m_rve, n, ke);

		// create the element's residual
		fe.assign(ndof, 0);

		// calculate the element's residual
		bd.ElementInternalForce(e, fe);

		// loop over the element's nodes
		for (i=0; i<ne; ++i)
		{
			FENode& ni = m.Node(e.m_node[i]);
			for (j=0; j<ne; ++j)
			{
				FENode& nj = m.Node(e.m_node[j]);
				if ((ni.m_ID[DOF_X] < 0) && (nj.m_ID[DOF_X] < 0))
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
				}
			}
/*
			if (ni.m_ID[DOF_X] < 0)
			{
				vec3d ri = ni.m_r0;

				double Fi[3] = {fe[3*i], fe[3*i+1], fe[3*i+2] };
				double Ri[3] = { ri.x, ri.y, ri.z };
				double I[3][3] = {{1,0,0},{0,1,0},{0,0,1}};

				D[0][0] += Fi[0]*Ri[0];
				D[1][1] += Fi[1]*Ri[1];
				D[2][2] += Fi[2]*Ri[2];

//				D[0][1] += 0;
//				D[0][2] += 0;
//				D[1][2] += 0;

				D[0][3] += 0.5*(Fi[0]*I[0][0]*Ri[1] + Fi[0]*I[0][1]*Ri[0]);
				D[0][4] += 0.5*(Fi[0]*I[0][1]*Ri[2] + Fi[0]*I[0][2]*Ri[1]);
				D[0][5] += 0.5*(Fi[0]*I[0][0]*Ri[2] + Fi[0]*I[0][2]*Ri[0]);

				D[1][3] += 0.5*(Fi[1]*I[1][0]*Ri[1] + Fi[1]*I[1][1]*Ri[0]);
				D[1][4] += 0.5*(Fi[1]*I[1][1]*Ri[2] + Fi[1]*I[1][2]*Ri[1]);
				D[1][5] += 0.5*(Fi[1]*I[1][0]*Ri[2] + Fi[1]*I[1][2]*Ri[0]);

				D[2][3] += 0.5*(Fi[2]*I[2][0]*Ri[1] + Fi[2]*I[2][1]*Ri[0]);
				D[2][4] += 0.5*(Fi[2]*I[2][1]*Ri[2] + Fi[2]*I[2][2]*Ri[1]);
				D[2][5] += 0.5*(Fi[2]*I[2][0]*Ri[2] + Fi[2]*I[2][2]*Ri[0]);

				D[3][3] += 0.25*(Fi[0]*I[1][0]*Ri[1] + Fi[1]*I[0][0]*Ri[1] + Fi[0]*I[1][1]*Ri[0] + Fi[1]*I[0][1]*Ri[0]);
				D[3][4] += 0.25*(Fi[0]*I[1][1]*Ri[2] + Fi[1]*I[0][1]*Ri[2] + Fi[0]*I[1][2]*Ri[1] + Fi[1]*I[0][2]*Ri[1]);
				D[3][5] += 0.25*(Fi[0]*I[1][0]*Ri[2] + Fi[1]*I[0][0]*Ri[2] + Fi[0]*I[1][2]*Ri[0] + Fi[1]*I[0][2]*Ri[0]);

				D[4][4] += 0.25*(Fi[1]*I[2][1]*Ri[2] + Fi[2]*I[1][1]*Ri[2] + Fi[1]*I[2][2]*Ri[1] + Fi[2]*I[1][2]*Ri[1]);
				D[4][5] += 0.25*(Fi[1]*I[2][0]*Ri[2] + Fi[2]*I[1][0]*Ri[2] + Fi[1]*I[2][2]*Ri[0] + Fi[2]*I[1][2]*Ri[0]);

				D[5][5] += 0.25*(Fi[0]*I[2][0]*Ri[2] + Fi[2]*I[0][0]*Ri[2] + Fi[0]*I[2][2]*Ri[0] + Fi[2]*I[0][2]*Ri[0]);
			}
*/
		}
	}

	// divide by volume
	double Vi = 1.0/(pt.m_J * m_V0);
	D[0][0] *= Vi; D[0][1] *= Vi; D[0][2] *= Vi; D[0][3] *= Vi; D[0][4] *= Vi; D[0][5] *= Vi;
	D[1][1] *= Vi; D[1][2] *= Vi; D[1][3] *= Vi; D[1][4] *= Vi; D[1][5] *= Vi;
	D[2][2] *= Vi; D[2][3] *= Vi; D[2][4] *= Vi; D[2][5] *= Vi;
	D[3][3] *= Vi; D[3][4] *= Vi; D[3][5] *= Vi;
	D[4][4] *= Vi; D[4][5] *= Vi;
	D[5][5] *= Vi;
/*
	// add the stress contribution
	D[0][0] -= s.xx()*Fi[0][0];
	D[1][1] -= s.yy()*Fi[1][1];
	D[2][2] -= s.zz()*Fi[2][2];

	D[0][1] -= s.xx()*Fi[1][1];
	D[0][2] -= s.xx()*Fi[2][2];
	D[1][2] -= s.yy()*Fi[2][2];

	D[0][3] -= 0.5*s.xx()*(Fi[0][1] + Fi[1][0]);
	D[0][4] -= 0.5*s.xx()*(Fi[1][2] + Fi[2][1]);
	D[0][5] -= 0.5*s.xx()*(Fi[0][2] + Fi[2][0]);

	D[1][3] -= 0.5*s.yy()*(Fi[0][1] + Fi[1][0]);
	D[1][4] -= 0.5*s.yy()*(Fi[1][2] + Fi[2][1]);
	D[1][5] -= 0.5*s.yy()*(Fi[0][2] + Fi[2][0]);

	D[2][3] -= 0.5*s.zz()*(Fi[0][1] + Fi[1][0]);
	D[2][4] -= 0.5*s.zz()*(Fi[1][2] + Fi[2][1]);
	D[2][5] -= 0.5*s.zz()*(Fi[0][2] + Fi[2][0]);

	D[3][3] -= 0.5*s.xy()*(Fi[0][1] + Fi[1][0]);
	D[3][4] -= 0.5*s.xy()*(Fi[1][2] + Fi[2][1]);
	D[3][5] -= 0.5*s.xy()*(Fi[2][0] + Fi[0][2]);

	D[4][4] -= 0.5*s.yz()*(Fi[1][2] + Fi[2][1]);
	D[4][5] -= 0.5*s.yz()*(Fi[0][2] + Fi[2][0]);

	D[5][5] -= 0.5*s.xz()*(Fi[0][2] + Fi[2][0]);
*/

	return tens4ds(D);
}
