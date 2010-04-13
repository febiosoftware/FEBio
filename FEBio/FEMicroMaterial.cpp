#include "stdafx.h"
#include "FEMicroMaterial.h"
#include "FEElemElemList.h"
#include "log.h"
#include "console.h"
#include "FESolidSolver.h"

// register the material with the framework
REGISTER_MATERIAL(FEMicroMaterial, "micro-material");

// define the material parameters
BEGIN_PARAMETER_LIST(FEMicroMaterial, FESolidMaterial)
	ADD_PARAMETER(m_szrve, FE_PARAM_STRING, "RVE");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEMicroMaterial::FEMicroMaterial(void)
{
	m_szrve[0] = 0;
}

//-----------------------------------------------------------------------------
FEMicroMaterial::~FEMicroMaterial(void)
{
}

//-----------------------------------------------------------------------------
void FEMicroMaterial::Init()
{
	// try to load the RVE model
	if (m_rve.Input(m_szrve) == false)
	{
		throw MaterialError("An error occured trying to read the RVE model from file %s.", m_szrve);
	}

	// make sure the RVE problem doesn't output anything to a plot file
	m_rve.m_pStep->SetPlotLevel(FE_PLOT_NEVER);

	// make sure we are using the same linear solver as the parent FEM
	m_rve.m_nsolver = SKYLINE_SOLVER;

	// create the DC's for this RVE
	PrepRVE();
}

//-----------------------------------------------------------------------------
void FEMicroMaterial::PrepRVE()
{
	// first we need to find all the boundary nodes
	FEMesh& m = m_rve.m_mesh;
	int N = m.Nodes();
	vector<int> tag(N); tag.zero();

	// create the element-element list
	FEElemElemList EEL;
	EEL.Create(&m);

	// use the E-E list to tag all exterior nodes
	int fn[4], nf, M = 0;
	FEElasticSolidDomain& bd = dynamic_cast<FEElasticSolidDomain&>(m.Domain(0));
	for (int i=0; i<bd.Elements(); ++i, ++M)
	{
		FESolidElement& el = bd.Element(i);
		nf = m.Faces(el);
		for (int j=0; j<nf; ++j)
		{
			if (EEL.Neighbor(M, j) == 0)
			{
				int nn = m.GetFace(el, j, fn);

				// mark all nodes
				tag[ fn[0] ] = 1;
				tag[ fn[1] ] = 1;
				tag[ fn[2] ] = 1;
				if (nn == 4) tag[ fn[3] ] = 1;
			}
		}
	}

	// count the nr of exterior nodes
	int NN = 0, i;
	for (i=0; i<N; ++i) if (tag[i] == 1) ++NN;

	assert(NN > 0);

	// create the DC's
	m_rve.m_DC.setsize(NN*3);
	NN = 0;
	for (i=0; i<N; ++i)
		if (tag[i] == 1)
		{
			for (int j=0; j<3; ++j, ++NN)
			{
				FENodalDisplacement& dc = m_rve.m_DC[NN];
				dc.bc = j;
				dc.lc = 0;	// we use the zeroth loadcurve
				dc.node = i;
				dc.s = 0;
			}
		}

	// the logfile is a shared resource between the master FEM and the RVE
	// in order not to corrupt the logfile we don't print anything for
	// the RVE problem.
	Logfile& log = GetLogfile();
	Logfile::MODE nmode = log.GetMode();
	GetLogfile().SetMode(Logfile::NEVER);

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
		bd.UnpackElement(el);
		nint = el.GaussPoints();
		w = el.GaussWeights();
		ve = 0;
		for (n=0; n<nint; ++n)
		{
			FEElasticMaterialPoint& pt = *el.m_State[n]->ExtractData<FEElasticMaterialPoint>();
			J = el.detJt(n);

			ve += J*w[n];
		}
		m_V0 += ve;
	}

	// reset the logfile mode
	log.SetMode(nmode);
}

//-----------------------------------------------------------------------------
mat3ds FEMicroMaterial::Stress(FEMaterialPoint &mp)
{
	// get the deformation gradient
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	mat3d F = pt.F;

	// the logfile is a shared resource between the master FEM and the RVE
	// in order not to corrupt the logfile we don't print anything for
	// the RVE problem.
	Logfile& log = GetLogfile();
	Logfile::MODE nmode = log.GetMode();
	GetLogfile().SetMode(Logfile::NEVER);

	// reset the RVE
	m_rve.Reset();

	// get the mesh
	FEMesh& m = m_rve.m_mesh;

	// assign new DC's for the boundary nodes
	int N = m_rve.m_DC.size()/3, i;
	for (i=0; i<N; ++i)
	{
		FENodalDisplacement& dx = m_rve.m_DC[3*i  ];
		FENodalDisplacement& dy = m_rve.m_DC[3*i+1];
		FENodalDisplacement& dz = m_rve.m_DC[3*i+2];

		FENode& node = m.Node(dx.node);

		vec3d r0 = node.m_r0;
		vec3d r1 = F*r0;

		dx.s = r1.x - r0.x;
		dy.s = r1.y - r0.y;
		dz.s = r1.z - r0.z;
	}

	// turn the console off
	Console::GetHandle()->Deactivate();

	// solve the RVE
	m_rve.Solve();

	// reset the logfile mode
	log.SetMode(nmode);

	// reactivate the console
	Console::GetHandle()->Activate();

	// calculate the averaged stress
	return AveragedStress(pt);
}

//-----------------------------------------------------------------------------

mat3ds FEMicroMaterial::AveragedStress(FEMaterialPoint& mp)
{
	// get the deformation gradient
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	mat3d F = pt.F;
	double J = pt.J;

	// get the mesh
	FEMesh& m = m_rve.m_mesh;
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
			J = el.detJt(n);

			ve += J*w[n];
			s += pt.s*(J*w[n]);
		}
		V += ve;
	}
	s /= V;
*/
	// get the reaction force vector from the solid solver
	FESolidSolver* ps = dynamic_cast<FESolidSolver*>(m_rve.m_pStep->m_psolver);
	assert(ps);
	vector<double>& R = ps->m_Fr;
	mat3ds s(0);
	for (int i=0; i<(int) m_rve.m_DC.size()/3; ++i)
	{
		FENodalDisplacement& dc = m_rve.m_DC[3*i];
		FENode& n = m.Node(dc.node);
		vec3d f;
		f.x = R[-n.m_ID[0]-2];
		f.y = R[-n.m_ID[1]-2];
		f.z = R[-n.m_ID[2]-2];
		s += (f & n.m_rt).sym();
	}
	s /= (J / m_V0);

	return s;
}

//-----------------------------------------------------------------------------
tens4ds FEMicroMaterial::Tangent(FEMaterialPoint &mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// get the mesh
	FEMesh& m = m_rve.m_mesh;

	// get the solver
	FESolidSolver* ps = dynamic_cast<FESolidSolver*>(m_rve.m_pStep->m_psolver);

	// the element's stiffness matrix
	matrix ke;

	// elasticity tensor
	double D[6][6] = {0};

	// get deformation gradient and its inverse
	mat3d F = pt.F;
	mat3d Fi = F.inverse();

	// get the stress
	mat3ds s = pt.s;

	// calculate the stiffness matrix
	FEElasticSolidDomain& bd = dynamic_cast<FEElasticSolidDomain&>(m.Domain(0));
	int NS = bd.Elements(), i, j;
	for (int n=0; n<NS; ++n)
	{
		FESolidElement& e = bd.Element(n);
		bd.UnpackElement(e);

		// create the element's stiffness matrix
		int ne = e.Nodes();
		int ndof = 3*ne;
		ke.Create(ndof, ndof);
		ke.zero();

		// calculate the element's stiffness matrix
		bd.ElementStiffness(m_rve, e, ke);

		// loop over the element's nodes
		for (i=0; i<ne; ++i)
		{
			FENode& ni = m.Node(e.m_node[i]);
			for (j=0; j<ne; ++j)
			{
				FENode& nj = m.Node(e.m_node[j]);
				if ((ni.m_ID[0] < 0) && (nj.m_ID[0] < 0))
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
				}
			}
		}
	}

	// divide by volume
	double Vi = pt.J / m_V0;
	D[0][0] *= Vi; D[0][1] *= Vi; D[0][2] *= Vi; D[0][3] *= Vi; D[0][4] *= Vi; D[0][5] *= Vi;
	D[1][1] *= Vi; D[1][2] *= Vi; D[1][3] *= Vi; D[1][4] *= Vi; D[1][5] *= Vi;
	D[2][2] *= Vi; D[2][3] *= Vi; D[2][4] *= Vi; D[2][5] *= Vi;
	D[3][3] *= Vi; D[3][4] *= Vi; D[3][5] *= Vi;
	D[4][4] *= Vi; D[4][5] *= Vi;
	D[5][5] *= Vi;

	return tens4ds(D);
}
