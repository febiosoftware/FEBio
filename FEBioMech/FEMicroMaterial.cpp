#include "stdafx.h"
#include "FEMicroMaterial.h"
#include "FECore/FEElemElemList.h"
#include "FECore/log.h"
#include "FESolidSolver.h"
#include "FEElasticSolidDomain.h"
#include "FECore/FEAnalysis.h"
#include "FEBioXML/FEBioImport.h"
#include "FEBioPlot/FEBioPlotFile.h"

//-----------------------------------------------------------------------------
FEMicroMaterialPoint::FEMicroMaterialPoint(FEMaterialPoint* mp) : FEMaterialPoint(mp)
{
}

//-----------------------------------------------------------------------------
//! Initialize material point data
void FEMicroMaterialPoint::Init(bool bflag)
{
}

//-----------------------------------------------------------------------------
//! create a shallow copy
FEMaterialPoint* FEMicroMaterialPoint::Copy()
{
	FEMicroMaterialPoint* pt = new FEMicroMaterialPoint(m_pt?m_pt->Copy():0);
	pt->m_Ka = m_Ka;
	return pt;
}

//-----------------------------------------------------------------------------
//! serialize material point data
void FEMicroMaterialPoint::Serialize(DumpFile& ar)
{
	if (ar.IsSaving())
	{
		ar << m_Ka;
	}
	else
	{
		ar >> m_Ka;
	}
}

//-----------------------------------------------------------------------------
//! stream material point data
void FEMicroMaterialPoint::ShallowCopy(DumpStream& dmp, bool bsave)
{
	if (bsave)
	{
		dmp << m_Ka;
	}
	else
	{
		dmp >> m_Ka;
	}
}

//=============================================================================

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_PARAMETER_LIST(FEMicroMaterial, FEElasticMaterial)
	ADD_PARAMETER(m_szrve    , FE_PARAM_STRING, "RVE"     );
	ADD_PARAMETER(m_szbc     , FE_PARAM_STRING, "bc_set"  );
	ADD_PARAMETER(m_bperiodic, FE_PARAM_BOOL  , "periodic");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEMicroMaterial::FEMicroMaterial(FEModel* pfem) : FEElasticMaterial(pfem)
{
	m_brve = false;

	// initialize parameters
	m_szrve[0] = 0;
	m_szbc[0] = 0;
	m_bperiodic = false;
}

//-----------------------------------------------------------------------------
FEMicroMaterial::~FEMicroMaterial(void)
{
}

//-----------------------------------------------------------------------------
FEMaterialPoint* FEMicroMaterial::CreateMaterialPointData()
{
	return new FEMicroMaterialPoint(new FEElasticMaterialPoint);
}

//-----------------------------------------------------------------------------
void FEMicroMaterial::Init()
{
	// try to load the RVE model
	if (m_brve == false)
	{
		// load the RVE model
		FEFEBioImport fim;
		if (fim.Load(m_rve, m_szrve) == false)
		{
			throw MaterialError("An error occured trying to read the RVE model from file %s.", m_szrve);
		}

		// set the pardiso solver as default
		m_rve.m_nsolver = PARDISO_SOLVER;

		// make sure the RVE problem doesn't output anything to a plot file
		m_rve.GetCurrentStep()->SetPlotLevel(FE_PLOT_NEVER);

		// create the BC's for this RVE
		if (PrepRVE() == false) throw MaterialError("An error occurred preparing RVE model");

		// mark that we read and processed the RVE successfully
		m_brve = true;
	}
}

//-----------------------------------------------------------------------------
bool FEMicroMaterial::PrepRVE()
{
	// find all boundar nodes
	FindBoundaryNodes();

	if (m_bperiodic == false)
	{
		// prep displacement BC's
		if (PrepDisplacementBC() == false) return false;
	}
	else
	{
		// prep periodic BC's
		if (PrepPeriodicBC() == false) return false;
	}

	// the logfile is a shared resource between the master FEM and the RVE
	// in order not to corrupt the logfile we don't print anything for
	// the RVE problem.
	Logfile::MODE nmode = felog.GetMode();
	felog.SetMode(Logfile::NEVER);

	// initialize RVE
	if (m_rve.Init() == false) return false;

	// calculate intial RVE volume
	m_V0 = 0;
	double ve;
	int nint;
	double* w, J;
	FEMesh& m = m_rve.GetMesh();
	for (int k=0; k<m.Domains(); ++k)
	{
		FESolidDomain& dom = static_cast<FESolidDomain&>(m.Domain(k));
		for (int i=0; i<dom.Elements(); ++i)
		{
			FESolidElement& el = dom.Element(i);
			nint = el.GaussPoints();
			w = el.GaussWeights();
			ve = 0;
			for (int n=0; n<nint; ++n)
			{
				FEElasticMaterialPoint& pt = *el.GetMaterialPoint(n)->ExtractData<FEElasticMaterialPoint>();
				J = dom.detJt(el, n);

				ve += J*w[n];
			}
			m_V0 += ve;
		}
	}

	// reset the logfile mode
	felog.SetMode(nmode);

	return true;
}

//-----------------------------------------------------------------------------
void FEMicroMaterial::FindBoundaryNodes()
{
	// first we need to find all the boundary nodes
	FEMesh& m = m_rve.GetMesh();
	int N = m.Nodes();
	m_BN.assign(N, 0);

	// create the element-element list
	FEElemElemList EEL;
	EEL.Create(&m);

	// use the E-E list to tag all exterior nodes
	int fn[FEElement::MAX_NODES], nf, M = 0;
	for (int k=0; k<m.Domains(); ++k)
	{
		FEDomain& dom = m.Domain(k);
		for (int i=0; i<dom.Elements(); ++i, ++M)
		{
			FEElement& el = dom.ElementRef(i);
			nf = m.Faces(el);
			for (int j=0; j<nf; ++j)
			{
				if (EEL.Neighbor(M, j) == 0)
				{
					// mark all nodes
					int nn = m.GetFace(el, j, fn);
					for (int k=0; k<nn; ++k) m_BN[fn[k]] = 1;
				}
			}
		}
	}

}

//-----------------------------------------------------------------------------
bool FEMicroMaterial::PrepDisplacementBC()
{
	FEMesh& m = m_rve.GetMesh();
	int N = m.Nodes();

	// count the nr of exterior nodes
	int NN = 0, i;
	for (i=0; i<N; ++i) if (m_BN[i] == 1) ++NN;

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
		if (m_BN[i] == 1)
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

	return true;
}

//-----------------------------------------------------------------------------
bool FEMicroMaterial::PrepPeriodicBC()
{
	// get the RVE mesh
	FEMesh& m = m_rve.GetMesh();

	// create a load curve
	FELoadCurve* plc = new FELoadCurve;
	plc->SetInterpolation(FELoadCurve::LINEAR);
	plc->Add(0.0, 0.0);
	plc->Add(1.0, 1.0);
	m_rve.AddLoadCurve(plc);
	int NLC = m_rve.LoadCurves() - 1;

	// find the node set that defines the corner nodes
	FENodeSet* pset = m.FindNodeSet(m_szbc);
	if (pset == 0) return false;

	// create the DC's
	m_rve.ClearBCs();
	int N = pset->size();
	for (int i=0; i<N; ++i)
		for (int j=0; j<3; ++j)
		{
			FEPrescribedBC* pdc = new FEPrescribedBC(&m_rve);
			pdc->bc = j;
			pdc->lc = NLC;
			pdc->node = (*pset)[i];
			pdc->s = 0;
			m_rve.AddPrescribedBC(pdc);
		}

	return true;
}

//-----------------------------------------------------------------------------
//! Assign the prescribed displacement to the boundary nodes.
void FEMicroMaterial::UpdateBC(FEModel& rve, mat3d& F)
{
	// get the mesh
	FEMesh& m = rve.GetMesh();

	// assign new DC's for the boundary nodes
	int N = rve.PrescribedBCs()/3, i;
	for (i=0; i<N; ++i)
	{
		FEPrescribedBC& dx = *rve.PrescribedBC(3*i  );
		FEPrescribedBC& dy = *rve.PrescribedBC(3*i+1);
		FEPrescribedBC& dz = *rve.PrescribedBC(3*i+2);

		FENode& node = m.Node(dx.node);

		vec3d r0 = node.m_r0;
		vec3d r1 = F*r0;

		dx.s = r1.x - r0.x;
		dy.s = r1.y - r0.y;
		dz.s = r1.z - r0.z;
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
			FEPeriodicBoundary* pc = dynamic_cast<FEPeriodicBoundary*>(rve.SurfacePairInteraction(i));
			assert(pc);

			// get the position of the first node
			vec3d r0 = pc->m_ss.Node(0).m_r0;

			// calculate the position of the projection
			FESurfaceElement* pm = pc->m_ss.m_pme[0]; assert(pm);
			for (int j=0; j<pm->Nodes(); ++j) r[j] = m.Node(pm->m_node[j]).m_r0;
			vec2d q = pc->m_ss.m_rs[0];
			vec3d r1 = pm->eval(r, q[0], q[1]);

			// calculate the offset distance
			vec3d u0 = r1 - r0;

			// apply deformation
			vec3d u1 = U*u0;

			// set this as the scale parameter for the offset
			FEParam* pp = pc->GetParameterList().Find("offset");
			assert(pp);
			pp->m_vscl = u1;
			pp->m_nlc = 0;
		}
	}
}

//-----------------------------------------------------------------------------
mat3ds FEMicroMaterial::Stress(FEMaterialPoint &mp)
{
	// get the deformation gradient
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	FEMicroMaterialPoint& mmpt = *mp.ExtractData<FEMicroMaterialPoint>();
	mat3d F = pt.m_F;

	// Create a local copy of the rve
	FEModel rve;
	rve.CopyFrom(m_rve);
	rve.GetStep(0)->SetPrintLevel(FE_PRINT_NEVER);

	// initialize
	if (rve.Init() == false) throw FEMultiScaleException();

	// apply the BC's
	UpdateBC(rve, F);

	// solve the RVE
	bool bret = rve.Solve();
/*
	// set the plot file
	FEBioPlotFile* pplt = new FEBioPlotFile(rve);
	vector<int> item;
	pplt->AddVariable("displacement", item);
	pplt->AddVariable("stress", item);
	pplt->Open(rve, "rve.xplt");
	pplt->Write(rve);
	pplt->Close();
*/
	// make sure it converged
	if (bret == false) throw FEMultiScaleException();

	// calculate the averaged stress
	mat3ds sa = AveragedStress(rve, mp);

	// calculate the averaged stiffness
	mmpt.m_Ka = AveragedStiffness(rve, mp);

	return sa;
}

//-----------------------------------------------------------------------------
// The stiffness is evaluated at the same time the stress is evaluated so we 
// can just return it here. Note that this assumes that the stress function 
// is always called prior to the tangent function.
tens4ds FEMicroMaterial::Tangent(FEMaterialPoint &mp)
{
	FEMicroMaterialPoint& mmpt = *mp.ExtractData<FEMicroMaterialPoint>();
	return mmpt.m_Ka;
}

//-----------------------------------------------------------------------------
//! Calculate the average stress from the RVE solution.
mat3ds FEMicroMaterial::AveragedStress(FEModel& rve, FEMaterialPoint &mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	mat3d F = pt.m_F;
	double J = pt.m_J;

	// get the RVE mesh
	FEMesh& m = rve.GetMesh();

	mat3d T; T.zero();

	// for periodic BC's we take the reaction forces directly from the periodic constraints
	if (m_bperiodic)
	{
		// get the reaction for from the periodic constraints
		for (int i=0; i<3; ++i)
		{
			FEPeriodicBoundary* pbc = dynamic_cast<FEPeriodicBoundary*>(rve.SurfacePairInteraction(i));
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
				T += (f & node.m_rt)*2.0;
			}
		}
	}

	// get the reaction force vector from the solid solver
	// (We also need to do this for the periodic BC, since at the prescribed nodes,
	// the contact forces will be zero). 
	FEAnalysis* pstep = rve.GetCurrentStep();
	FESolidSolver* ps = dynamic_cast<FESolidSolver*>(pstep->m_psolver);
	assert(ps);
	vector<double>& R = ps->m_Fr;
	int nbc = rve.PrescribedBCs();
	for (int i=0; i<nbc/3; ++i)
	{
		FEPrescribedBC& dc = *rve.PrescribedBC(3*i);
		FENode& n = m.Node(dc.node);
		vec3d f;
		f.x = R[-n.m_ID[DOF_X]-2];
		f.y = R[-n.m_ID[DOF_Y]-2];
		f.z = R[-n.m_ID[DOF_Z]-2];
		T += f & n.m_rt;
	}

	return T.sym() / (J*m_V0);
}

//-----------------------------------------------------------------------------
//! Calculate the average stiffness from the RVE solution.
tens4ds FEMicroMaterial::AveragedStiffness(FEModel& rve, FEMaterialPoint &mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// get the mesh
	FEMesh& m = rve.GetMesh();

	// get the solver
	FEAnalysis* pstep = rve.GetCurrentStep();
	FESolidSolver* ps = dynamic_cast<FESolidSolver*>(pstep->m_psolver);

	// the element's stiffness matrix
	matrix ke;

	// element's residual
	vector<double> fe;

	// get deformation gradient and its inverse
	mat3d F = pt.m_F;
	mat3d Fi = F.inverse();

	// get the stress
	mat3ds s = pt.m_s;

	// calculate the center point
	vec3d rc(0,0,0);
	for (int k=0; k<m.Nodes(); ++k) rc += m.Node(k).m_rt;
	rc /= (double) m.Nodes();

	// LTE - Calculate the initial center point
	vec3d rc0(0,0,0);
	for (int k=0; k<m.Nodes(); ++k) rc0 += m.Node(k).m_r0;
	rc0 /= (double) m.Nodes();

	// LTE - 
	bool old_flag = true;
	tens4ds c(0.);

	if (old_flag)
	{
		// LTE - elasticity tensor
		double D[6][6] = {0};
		
		// calculate the stiffness matrix and residual
		for (int k=0; k<m.Domains(); ++k)
		{
			FEElasticSolidDomain& bd = dynamic_cast<FEElasticSolidDomain&>(m.Domain(k));
			int NS = bd.Elements();
			for (int n=0; n<NS; ++n)
			{
				FESolidElement& e = bd.Element(n);

				// create the element's stiffness matrix
				int ne = e.Nodes();
				int ndof = 3*ne;
				ke.resize(ndof, ndof);
				ke.zero();

				// calculate the element's stiffness matrix
				bd.ElementStiffness(rve, n, ke);

				// create the element's residual
				fe.assign(ndof, 0);

				// calculate the element's residual
				bd.ElementInternalForce(e, fe);

				// loop over the element's nodes
				for (int i=0; i<ne; ++i)
				{
					FENode& ni = m.Node(e.m_node[i]);
					for (int j=0; j<ne; ++j)
					{
						FENode& nj = m.Node(e.m_node[j]);
						if ((m_BN[e.m_node[i]] == 1) && (m_BN[e.m_node[j]] == 1))
						{
							// both nodes are boundary nodes
							// so grab the element's submatrix
							double K[3][3];
							K[0][0] = ke[3*i  ][3*j  ]; K[0][1] = ke[3*i  ][3*j+1]; K[0][2] = ke[3*i  ][3*j+2];
							K[1][0] = ke[3*i+1][3*j  ]; K[1][1] = ke[3*i+1][3*j+1]; K[1][2] = ke[3*i+1][3*j+2];
							K[2][0] = ke[3*i+2][3*j  ]; K[2][1] = ke[3*i+2][3*j+1]; K[2][2] = ke[3*i+2][3*j+2];

							// get the nodal positions
							vec3d ri = ni.m_rt - rc;
							vec3d rj = nj.m_rt - rc;

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
	
		c = tens4ds(D);
	}	
	else
	{
		double A[9][9] = {0.};
		
		// calculate the stiffness matrix and residual
		for (int k=0; k<m.Domains(); ++k)
		{
			FEElasticSolidDomain& bd = dynamic_cast<FEElasticSolidDomain&>(m.Domain(k));
			int NS = bd.Elements();
			for (int n=0; n<NS; ++n)
			{
				FESolidElement& e = bd.Element(n);

				// create the element's stiffness matrix
				int ne = e.Nodes();
				int ndof = 3*ne;
				ke.resize(ndof, ndof);
				ke.zero();

				// calculate the element's stiffness matrix
				bd.ElementStiffness(rve, n, ke);

				// create the element's residual
				fe.assign(ndof, 0);

				// calculate the element's residual
				bd.ElementInternalForce(e, fe);

				// loop over the element's nodes
				for (int i=0; i<ne; ++i)
				{
					FENode& ni = m.Node(e.m_node[i]);
					for (int j=0; j<ne; ++j)
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

							// LTE - get the initial nodal positions
							vec3d ri0 = ni.m_r0 - rc0;
							vec3d rj0 = nj.m_r0 - rc0;
							
							double Ri[3] = { ri0.x, ri0.y, ri0.z };
							double Rj[3] = { rj0.x, rj0.y, rj0.z };

							// LTE - Create the first elasticity tensor
							// Calculate based on equation from Kouznetsova thesis 2.39
							// Aabcd = (1/V0)*XbKacXd
							// Major symmetrize in order to ensure that this tensor has major symmetry (we'll handle minor symmetry later)	
							// (1/2)*(Aabcd + Acdab) = (1/V0)*(1/2)*(XbKacXd + XdKcaXb)  
							A[0][0] += Ri[0]*K[0][0]*Rj[0]; 
							A[0][1] += 0.5*(Ri[0]*K[0][1]*Rj[1] + Ri[1]*K[1][0]*Rj[0]);
							A[0][2] += 0.5*(Ri[0]*K[0][2]*Rj[2] + Ri[2]*K[2][0]*Rj[0]);
							A[0][3] += 0.5*(Ri[0]*K[0][0]*Rj[1] + Ri[1]*K[0][0]*Rj[0]);
							A[0][4] += 0.5*(Ri[0]*K[0][1]*Rj[2] + Ri[2]*K[1][0]*Rj[0]);
							A[0][5] += 0.5*(Ri[0]*K[0][0]*Rj[2] + Ri[2]*K[0][0]*Rj[0]);
							A[0][6] += 0.5*(Ri[0]*K[0][1]*Rj[0] + Ri[0]*K[1][0]*Rj[0]);
							A[0][7] += 0.5*(Ri[0]*K[0][2]*Rj[1] + Ri[1]*K[2][0]*Rj[0]);
							A[0][8] += 0.5*(Ri[0]*K[0][2]*Rj[0] + Ri[0]*K[2][0]*Rj[0]);
							
							A[1][1] += Ri[1]*K[1][1]*Rj[1]; 
							A[1][2] += 0.5*(Ri[1]*K[1][2]*Rj[2] + Ri[2]*K[2][1]*Rj[1]); 
							A[1][3] += 0.5*(Ri[1]*K[1][0]*Rj[1] + Ri[1]*K[0][1]*Rj[1]);
							A[0][4] += 0.5*(Ri[1]*K[1][1]*Rj[2] + Ri[2]*K[1][1]*Rj[1]);
							A[0][5] += 0.5*(Ri[1]*K[1][0]*Rj[2] + Ri[2]*K[0][1]*Rj[1]);
							A[0][6] += 0.5*(Ri[1]*K[1][1]*Rj[0] + Ri[0]*K[1][1]*Rj[1]);
							A[0][7] += 0.5*(Ri[1]*K[1][2]*Rj[1] + Ri[1]*K[2][1]*Rj[1]);
							A[0][8] += 0.5*(Ri[1]*K[1][2]*Rj[0] + Ri[0]*K[2][1]*Rj[1]);

							A[2][2] += Ri[2]*K[2][2]*Rj[2]; 
							A[2][3] += 0.5*(Ri[2]*K[2][0]*Rj[1] + Ri[1]*K[0][2]*Rj[2]);
							A[2][4] += 0.5*(Ri[2]*K[2][1]*Rj[2] + Ri[2]*K[1][2]*Rj[2]);
							A[2][5] += 0.5*(Ri[2]*K[2][0]*Rj[2] + Ri[2]*K[0][2]*Rj[2]);
							A[2][6] += 0.5*(Ri[2]*K[2][1]*Rj[0] + Ri[0]*K[1][2]*Rj[2]);
							A[2][7] += 0.5*(Ri[2]*K[2][2]*Rj[1] + Ri[1]*K[2][2]*Rj[2]);
							A[2][8] += 0.5*(Ri[2]*K[2][2]*Rj[0] + Ri[0]*K[2][2]*Rj[2]);

							A[3][3] += Ri[1]*K[0][0]*Rj[1];
							A[3][4] += 0.5*(Ri[1]*K[0][1]*Rj[2] + Ri[2]*K[1][0]*Rj[1]);
							A[3][5] += 0.5*(Ri[1]*K[0][0]*Rj[2] + Ri[2]*K[0][0]*Rj[1]);
							A[3][6] += 0.5*(Ri[1]*K[0][1]*Rj[0] + Ri[0]*K[1][0]*Rj[1]);
							A[3][7] += 0.5*(Ri[1]*K[0][2]*Rj[1] + Ri[1]*K[2][0]*Rj[1]);
							A[3][8] += 0.5*(Ri[1]*K[0][2]*Rj[0] + Ri[0]*K[2][0]*Rj[1]);

							A[4][4] += Ri[2]*K[1][1]*Rj[2];
							A[4][5] += 0.5*(Ri[2]*K[1][0]*Rj[2] + Ri[2]*K[0][1]*Rj[2]);
							A[4][6] += 0.5*(Ri[2]*K[1][1]*Rj[0] + Ri[0]*K[1][1]*Rj[2]);
							A[4][7] += 0.5*(Ri[2]*K[1][2]*Rj[1] + Ri[1]*K[2][1]*Rj[2]);
							A[4][8] += 0.5*(Ri[2]*K[1][2]*Rj[0] + Ri[0]*K[2][1]*Rj[2]);

							A[5][5] += Ri[2]*K[0][0]*Rj[2];
							A[5][6] += 0.5*(Ri[2]*K[0][1]*Rj[0] + Ri[0]*K[1][0]*Rj[2]);
							A[5][7] += 0.5*(Ri[2]*K[0][2]*Rj[1] + Ri[1]*K[2][0]*Rj[2]);
							A[5][8] += 0.5*(Ri[2]*K[0][2]*Rj[0] + Ri[0]*K[2][0]*Rj[2]);

							A[6][6] += Ri[0]*K[1][1]*Rj[0];
							A[6][7] += 0.5*(Ri[0]*K[1][2]*Rj[1] + Ri[1]*K[2][1]*Rj[0]);
							A[6][8] += 0.5*(Ri[0]*K[1][2]*Rj[0] + Ri[0]*K[2][1]*Rj[0]);

							A[7][7] += Ri[1]*K[2][2]*Rj[1];
							A[7][8] += 0.5*(Ri[1]*K[2][2]*Rj[0] + Ri[0]*K[2][2]*Rj[1]);

							A[8][8] += Ri[0]*K[2][2]*Rj[0];

						}
					}																										
				}
			}
		}

		// LTE - divide by initial volume
		tens4dms a = tens4dms(A);
		a /= m_V0;
				
		// LTE - Calculate the major symmetric spatial elasticity tensor based on the following equation:
		// CMSabcd = (1/J)*Fbi*Aaicj*Fjd - IacSbd 
		// ?Maybe move this to its own function?
		tens4dms cms(0.);

        cms.d[0] = a.d[0]*F[0][0]*F[0][0] + a.d[15]*F[0][0]*F[0][2] + a.d[18]*F[0][2]*F[1][0] + a.d[15]*F[0][0]*F[2][0] + a.d[18]*F[0][1]*F[2][0] + a.d[20]*F[0][2]*F[2][0] + a.d[6]*F[0][0]*F[0][1] + a.d[6]*F[0][0]*F[1][0] + a.d[9]*F[0][1]*F[1][0];
		
		cms.d[1] = a.d[24]*F[0][1]*F[0][1] + a.d[1]*F[0][0]*F[1][1] + a.d[21]*F[0][0]*F[0][1] + a.d[16]*F[0][2]*F[1][1] + a.d[26]*F[0][1]*F[0][2] + a.d[10]*F[0][0]*F[2][1] + a.d[13]*F[0][1]*F[2][1] + a.d[19]*F[0][2]*F[2][1] + a.d[7]*F[0][1]*F[1][1];
		cms.d[2] = a.d[2]*F[1][1]*F[1][1] + a.d[11]*F[1][1]*F[1][2] + a.d[22]*F[0][1]*F[1][1] + a.d[25]*F[0][1]*F[1][2] + a.d[27]*F[0][1]*F[1][0] + a.d[11]*F[1][1]*F[2][1] + a.d[22]*F[1][0]*F[1][1] + a.d[14]*F[1][2]*F[2][1] + a.d[25]*F[1][0]*F[2][1];

		cms.d[3] = a.d[41]*F[0][2]*F[0][2] + a.d[36]*F[0][0]*F[0][2] + a.d[28]*F[0][0]*F[1][2] + a.d[17]*F[0][2]*F[2][2] + a.d[39]*F[0][1]*F[0][2] + a.d[31]*F[0][1]*F[1][2] + a.d[33]*F[0][2]*F[1][2] + a.d[3]*F[0][0]*F[2][2] + a.d[8]*F[0][1]*F[2][2];
		cms.d[4] = a.d[32]*F[1][2]*F[1][2] + a.d[12]*F[1][2]*F[2][2] + a.d[37]*F[0][2]*F[1][1] + a.d[29]*F[1][1]*F[1][2] + a.d[40]*F[0][2]*F[1][2] + a.d[42]*F[0][2]*F[1][0] + a.d[23]*F[1][0]*F[2][2] + a.d[34]*F[1][0]*F[1][2] + a.d[4]*F[1][1]*F[2][2];
		cms.d[5] = a.d[5]*F[2][2]*F[2][2] + a.d[38]*F[0][2]*F[2][2] + a.d[30]*F[1][2]*F[2][2] + a.d[43]*F[0][2]*F[2][1] + a.d[44]*F[0][2]*F[2][0] + a.d[35]*F[1][2]*F[2][1] + a.d[30]*F[2][1]*F[2][2] + a.d[43]*F[1][2]*F[2][0] + a.d[38]*F[2][0]*F[2][2];
		
		cms.d[6] = a.d[6]*F[0][1]*F[0][1] + a.d[0]*F[0][0]*F[0][1] + a.d[15]*F[0][1]*F[0][2] + a.d[18]*F[0][2]*F[1][1] + a.d[15]*F[0][0]*F[2][1] + a.d[18]*F[0][1]*F[2][1] + a.d[20]*F[0][2]*F[2][1] + a.d[6]*F[0][0]*F[1][1] + a.d[9]*F[0][1]*F[1][1];
		cms.d[7] = a.d[7]*F[1][1]*F[1][1] + a.d[1]*F[0][1]*F[1][1] + a.d[10]*F[0][1]*F[1][2] + a.d[21]*F[0][1]*F[1][0] + a.d[13]*F[1][1]*F[1][2] + a.d[24]*F[1][0]*F[1][1] + a.d[16]*F[1][1]*F[2][1] + a.d[19]*F[1][2]*F[2][1] + a.d[26]*F[1][0]*F[2][1];
		cms.d[8] = a.d[33]*F[2][1]*F[2][1] + a.d[28]*F[0][1]*F[2][1] + a.d[3]*F[0][1]*F[2][2] + a.d[36]*F[0][1]*F[2][0] + a.d[17]*F[2][1]*F[2][2] + a.d[31]*F[1][1]*F[2][1] + a.d[39]*F[1][1]*F[2][0] + a.d[41]*F[2][0]*F[2][1] + a.d[8]*F[1][1]*F[2][2];
		cms.d[9] = a.d[9]*F[1][1]*F[1][1] + a.d[0]*F[0][1]*F[1][0] + a.d[15]*F[0][1]*F[1][2] + a.d[18]*F[1][1]*F[1][2] + a.d[15]*F[1][0]*F[2][1] + a.d[18]*F[1][1]*F[2][1] + a.d[20]*F[1][2]*F[2][1] + a.d[6]*F[0][1]*F[1][1] + a.d[6]*F[1][0]*F[1][1];

		cms.d[10] = a.d[26]*F[0][2]*F[0][2] + a.d[1]*F[0][0]*F[1][2] + a.d[21]*F[0][0]*F[0][2] + a.d[24]*F[0][1]*F[0][2] + a.d[16]*F[0][2]*F[1][2] + a.d[10]*F[0][0]*F[2][2] + a.d[13]*F[0][1]*F[2][2] + a.d[19]*F[0][2]*F[2][2] + a.d[7]*F[0][1]*F[1][2];
		cms.d[11] = a.d[11]*F[1][2]*F[1][2] + a.d[22]*F[0][2]*F[1][1] + a.d[25]*F[0][2]*F[1][2] + a.d[27]*F[0][2]*F[1][0] + a.d[2]*F[1][1]*F[1][2] + a.d[11]*F[1][1]*F[2][2] + a.d[22]*F[1][0]*F[1][2] + a.d[14]*F[1][2]*F[2][2] + a.d[25]*F[1][0]*F[2][2];
		cms.d[12] = a.d[12]*F[2][2]*F[2][2] + a.d[23]*F[0][2]*F[2][2] + a.d[34]*F[0][2]*F[2][1] + a.d[29]*F[1][2]*F[2][1] + a.d[42]*F[0][2]*F[2][0] + a.d[37]*F[1][2]*F[2][0] + a.d[4]*F[1][2]*F[2][2] + a.d[32]*F[2][1]*F[2][2] + a.d[40]*F[2][0]*F[2][2];
		cms.d[13] = a.d[16]*F[1][2]*F[1][2] + a.d[1]*F[1][0]*F[1][2] + a.d[21]*F[0][2]*F[1][0] + a.d[24]*F[0][2]*F[1][1] + a.d[26]*F[0][2]*F[1][2] + a.d[10]*F[1][0]*F[2][2] + a.d[13]*F[1][1]*F[2][2] + a.d[19]*F[1][2]*F[2][2] + a.d[7]*F[1][1]*F[1][2];
		cms.d[14] = a.d[14]*F[2][2]*F[2][2] + a.d[11]*F[1][2]*F[2][2] + a.d[22]*F[0][2]*F[2][1] + a.d[25]*F[0][2]*F[2][2] + a.d[27]*F[0][2]*F[2][0] + a.d[2]*F[1][2]*F[2][1] + a.d[11]*F[2][1]*F[2][2] + a.d[22]*F[1][2]*F[2][0] + a.d[25]*F[2][0]*F[2][2];

		cms.d[15] = a.d[15]*F[0][2]*F[0][2] + a.d[0]*F[0][0]*F[0][2] + a.d[18]*F[0][2]*F[1][2] + a.d[15]*F[0][0]*F[2][2] + a.d[18]*F[0][1]*F[2][2] + a.d[20]*F[0][2]*F[2][2] + a.d[6]*F[0][1]*F[0][2] + a.d[6]*F[0][0]*F[1][2] + a.d[9]*F[0][1]*F[1][2];
		cms.d[16] = a.d[13]*F[1][2]*F[1][2] + a.d[1]*F[0][2]*F[1][1] + a.d[10]*F[0][2]*F[1][2] + a.d[21]*F[0][2]*F[1][0] + a.d[24]*F[1][0]*F[1][2] + a.d[16]*F[1][1]*F[2][2] + a.d[19]*F[1][2]*F[2][2] + a.d[26]*F[1][0]*F[2][2] + a.d[7]*F[1][1]*F[1][2];
		cms.d[17] = a.d[17]*F[2][2]*F[2][2] + a.d[28]*F[0][2]*F[2][1] + a.d[3]*F[0][2]*F[2][2] + a.d[36]*F[0][2]*F[2][0] + a.d[31]*F[1][2]*F[2][1] + a.d[39]*F[1][2]*F[2][0] + a.d[33]*F[2][1]*F[2][2] + a.d[41]*F[2][0]*F[2][2] + a.d[8]*F[1][2]*F[2][2];
		cms.d[18] = a.d[18]*F[1][2]*F[1][2] + a.d[0]*F[0][2]*F[1][0] + a.d[15]*F[0][2]*F[1][2] + a.d[15]*F[1][0]*F[2][2] + a.d[18]*F[1][1]*F[2][2] + a.d[20]*F[1][2]*F[2][2] + a.d[6]*F[0][2]*F[1][1] + a.d[6]*F[1][0]*F[1][2] + a.d[9]*F[1][1]*F[1][2];
		cms.d[19] = a.d[19]*F[2][2]*F[2][2] + a.d[1]*F[0][2]*F[2][1] + a.d[10]*F[0][2]*F[2][2] + a.d[21]*F[0][2]*F[2][0] + a.d[13]*F[1][2]*F[2][2] + a.d[24]*F[1][2]*F[2][0] + a.d[16]*F[2][1]*F[2][2] + a.d[26]*F[2][0]*F[2][2] + a.d[7]*F[1][2]*F[2][1];
		cms.d[20] = a.d[20]*F[2][2]*F[2][2] + a.d[0]*F[0][2]*F[2][0] + a.d[15]*F[0][2]*F[2][2] + a.d[18]*F[1][2]*F[2][2] + a.d[15]*F[2][0]*F[2][2] + a.d[18]*F[2][1]*F[2][2] + a.d[6]*F[0][2]*F[2][1] + a.d[6]*F[1][2]*F[2][0] + a.d[9]*F[1][2]*F[2][1];

		cms.d[21] = a.d[21]*F[0][0]*F[0][0] + a.d[1]*F[0][0]*F[1][0] + a.d[24]*F[0][0]*F[0][1] + a.d[16]*F[0][2]*F[1][0] + a.d[26]*F[0][0]*F[0][2] + a.d[10]*F[0][0]*F[2][0] + a.d[13]*F[0][1]*F[2][0] + a.d[19]*F[0][2]*F[2][0] + a.d[7]*F[0][1]*F[1][0];
		cms.d[22] = a.d[22]*F[1][0]*F[1][0] + a.d[11]*F[1][0]*F[1][2] + a.d[22]*F[0][0]*F[1][1] + a.d[25]*F[0][0]*F[1][2] + a.d[27]*F[0][0]*F[1][0] + a.d[2]*F[1][0]*F[1][1] + a.d[11]*F[1][1]*F[2][0] + a.d[14]*F[1][2]*F[2][0] + a.d[25]*F[1][0]*F[2][0];
		cms.d[23] = a.d[40]*F[2][0]*F[2][0] + a.d[23]*F[0][0]*F[2][2] + a.d[12]*F[2][0]*F[2][2] + a.d[34]*F[0][0]*F[2][1] + a.d[29]*F[1][0]*F[2][1] + a.d[42]*F[0][0]*F[2][0] + a.d[37]*F[1][0]*F[2][0] + a.d[4]*F[1][0]*F[2][2] + a.d[32]*F[2][0]*F[2][1];
		cms.d[24] = a.d[1]*F[1][0]*F[1][0] + a.d[21]*F[0][0]*F[1][0] + a.d[24]*F[0][0]*F[1][1] + a.d[16]*F[1][0]*F[1][2] + a.d[26]*F[0][0]*F[1][2] + a.d[10]*F[1][0]*F[2][0] + a.d[13]*F[1][1]*F[2][0] + a.d[19]*F[1][2]*F[2][0] + a.d[7]*F[1][0]*F[1][1];
		cms.d[25] = a.d[25]*F[2][0]*F[2][0] + a.d[11]*F[1][0]*F[2][2] + a.d[22]*F[0][0]*F[2][1] + a.d[25]*F[0][0]*F[2][2] + a.d[27]*F[0][0]*F[2][0] + a.d[2]*F[1][0]*F[2][1] + a.d[11]*F[2][0]*F[2][1] + a.d[22]*F[1][0]*F[2][0] + a.d[14]*F[2][0]*F[2][2];
		cms.d[26] = a.d[10]*F[2][0]*F[2][0] + a.d[1]*F[1][0]*F[2][0] + a.d[21]*F[0][0]*F[2][0] + a.d[24]*F[0][0]*F[2][1] + a.d[16]*F[1][0]*F[2][2] + a.d[26]*F[0][0]*F[2][2] + a.d[13]*F[2][0]*F[2][1] + a.d[19]*F[2][0]*F[2][2] + a.d[7]*F[1][0]*F[2][1];
		cms.d[27] = a.d[27]*F[0][0]*F[0][0] + a.d[11]*F[0][2]*F[1][0] + a.d[22]*F[0][0]*F[0][1] + a.d[25]*F[0][0]*F[0][2] + a.d[2]*F[0][1]*F[1][0] + a.d[11]*F[0][1]*F[2][0] + a.d[22]*F[0][0]*F[1][0] + a.d[14]*F[0][2]*F[2][0] + a.d[25]*F[0][0]*F[2][0];

		cms.d[28] = a.d[39]*F[0][1]*F[0][1] + a.d[36]*F[0][0]*F[0][1] + a.d[28]*F[0][0]*F[1][1] + a.d[17]*F[0][2]*F[2][1] + a.d[31]*F[0][1]*F[1][1] + a.d[41]*F[0][1]*F[0][2] + a.d[33]*F[0][2]*F[1][1] + a.d[3]*F[0][0]*F[2][1] + a.d[8]*F[0][1]*F[2][1];
		cms.d[29] = a.d[29]*F[1][1]*F[1][1] + a.d[12]*F[1][2]*F[2][1] + a.d[37]*F[0][1]*F[1][1] + a.d[40]*F[0][1]*F[1][2] + a.d[42]*F[0][1]*F[1][0] + a.d[23]*F[1][0]*F[2][1] + a.d[32]*F[1][1]*F[1][2] + a.d[34]*F[1][0]*F[1][1] + a.d[4]*F[1][1]*F[2][1];
		cms.d[30] = a.d[30]*F[2][1]*F[2][1] + a.d[38]*F[0][1]*F[2][2] + a.d[30]*F[1][1]*F[2][2] + a.d[43]*F[0][1]*F[2][1] + a.d[44]*F[0][1]*F[2][0] + a.d[35]*F[1][1]*F[2][1] + a.d[43]*F[1][1]*F[2][0] + a.d[38]*F[2][0]*F[2][1] + a.d[5]*F[2][1]*F[2][2];
		cms.d[31] = a.d[31]*F[1][1]*F[1][1] + a.d[36]*F[0][1]*F[1][0] + a.d[28]*F[1][0]*F[1][1] + a.d[17]*F[1][2]*F[2][1] + a.d[39]*F[0][1]*F[1][1] + a.d[41]*F[0][1]*F[1][2] + a.d[33]*F[1][1]*F[1][2] + a.d[3]*F[1][0]*F[2][1] + a.d[8]*F[1][1]*F[2][1];
		cms.d[32] = a.d[4]*F[2][1]*F[2][1] + a.d[12]*F[2][1]*F[2][2] + a.d[37]*F[0][1]*F[2][1] + a.d[29]*F[1][1]*F[2][1] + a.d[40]*F[0][1]*F[2][2] + a.d[42]*F[0][1]*F[2][0] + a.d[23]*F[2][0]*F[2][1] + a.d[32]*F[1][1]*F[2][2] + a.d[34]*F[1][1]*F[2][0];
		cms.d[33] = a.d[8]*F[2][1]*F[2][1] + a.d[36]*F[0][1]*F[2][0] + a.d[28]*F[1][1]*F[2][0] + a.d[17]*F[2][1]*F[2][2] + a.d[39]*F[0][1]*F[2][1] + a.d[31]*F[1][1]*F[2][1] + a.d[41]*F[0][1]*F[2][2] + a.d[33]*F[1][1]*F[2][2] + a.d[3]*F[2][0]*F[2][1];
		cms.d[34] = a.d[37]*F[0][1]*F[0][1] + a.d[12]*F[0][2]*F[2][1] + a.d[29]*F[0][1]*F[1][1] + a.d[40]*F[0][1]*F[0][2] + a.d[42]*F[0][0]*F[0][1] + a.d[23]*F[0][0]*F[2][1] + a.d[32]*F[0][2]*F[1][1] + a.d[34]*F[0][0]*F[1][1] + a.d[4]*F[0][1]*F[2][1];
		cms.d[35] = a.d[35]*F[1][1]*F[1][1] + a.d[38]*F[0][1]*F[1][2] + a.d[30]*F[1][1]*F[1][2] + a.d[43]*F[0][1]*F[1][1] + a.d[44]*F[0][1]*F[1][0] + a.d[30]*F[1][1]*F[2][1] + a.d[43]*F[1][0]*F[1][1] + a.d[38]*F[1][0]*F[2][1] + a.d[5]*F[1][2]*F[2][1];

		cms.d[36] = a.d[36]*F[0][0]*F[0][0] + a.d[28]*F[0][0]*F[1][0] + a.d[17]*F[0][2]*F[2][0] + a.d[39]*F[0][0]*F[0][1] + a.d[31]*F[0][1]*F[1][0] + a.d[41]*F[0][0]*F[0][2] + a.d[33]*F[0][2]*F[1][0] + a.d[3]*F[0][0]*F[2][0] + a.d[8]*F[0][1]*F[2][0];
		cms.d[37] = a.d[34]*F[1][0]*F[1][0] + a.d[12]*F[1][2]*F[2][0] + a.d[37]*F[0][0]*F[1][1] + a.d[29]*F[1][0]*F[1][1] + a.d[40]*F[0][0]*F[1][2] + a.d[42]*F[0][0]*F[1][0] + a.d[23]*F[1][0]*F[2][0] + a.d[32]*F[1][0]*F[1][2] + a.d[4]*F[1][1]*F[2][0];
		cms.d[38] = a.d[38]*F[2][0]*F[2][0] + a.d[38]*F[0][0]*F[2][2] + a.d[30]*F[1][0]*F[2][2] + a.d[43]*F[0][0]*F[2][1] + a.d[44]*F[0][0]*F[2][0] + a.d[35]*F[1][0]*F[2][1] + a.d[30]*F[2][0]*F[2][1] + a.d[43]*F[1][0]*F[2][0] + a.d[5]*F[2][0]*F[2][2];
		cms.d[39] = a.d[28]*F[1][0]*F[1][0] + a.d[36]*F[0][0]*F[1][0] + a.d[17]*F[1][2]*F[2][0] + a.d[39]*F[0][0]*F[1][1] + a.d[31]*F[1][0]*F[1][1] + a.d[41]*F[0][0]*F[1][2] + a.d[33]*F[1][0]*F[1][2] + a.d[3]*F[1][0]*F[2][0] + a.d[8]*F[1][1]*F[2][0];
		cms.d[40] = a.d[23]*F[2][0]*F[2][0] + a.d[12]*F[2][0]*F[2][2] + a.d[37]*F[0][0]*F[2][1] + a.d[29]*F[1][0]*F[2][1] + a.d[40]*F[0][0]*F[2][2] + a.d[42]*F[0][0]*F[2][0] + a.d[32]*F[1][0]*F[2][2] + a.d[34]*F[1][0]*F[2][0] + a.d[4]*F[2][0]*F[2][1];
		cms.d[41] = a.d[3]*F[2][0]*F[2][0] + a.d[36]*F[0][0]*F[2][0] + a.d[28]*F[1][0]*F[2][0] + a.d[17]*F[2][0]*F[2][2] + a.d[39]*F[0][0]*F[2][1] + a.d[31]*F[1][0]*F[2][1] + a.d[41]*F[0][0]*F[2][2] + a.d[33]*F[1][0]*F[2][2] + a.d[8]*F[2][0]*F[2][1];
		cms.d[42] = a.d[42]*F[0][0]*F[0][0] + a.d[12]*F[0][2]*F[2][0] + a.d[37]*F[0][0]*F[0][1] + a.d[29]*F[0][1]*F[1][0] + a.d[40]*F[0][0]*F[0][2] + a.d[23]*F[0][0]*F[2][0] + a.d[32]*F[0][2]*F[1][0] + a.d[34]*F[0][0]*F[1][0] + a.d[4]*F[0][1]*F[2][0];
		cms.d[43] = a.d[43]*F[1][0]*F[1][0] + a.d[38]*F[0][0]*F[1][2] + a.d[30]*F[1][0]*F[1][2] + a.d[43]*F[0][0]*F[1][1] + a.d[44]*F[0][0]*F[1][0] + a.d[35]*F[1][0]*F[1][1] + a.d[30]*F[1][1]*F[2][0] + a.d[38]*F[1][0]*F[2][0] + a.d[5]*F[1][2]*F[2][0];
		cms.d[44] = a.d[44]*F[0][0]*F[0][0] + a.d[38]*F[0][0]*F[0][2] + a.d[30]*F[0][2]*F[1][0] + a.d[43]*F[0][0]*F[0][1] + a.d[35]*F[0][1]*F[1][0] + a.d[30]*F[0][1]*F[2][0] + a.d[43]*F[0][0]*F[1][0] + a.d[38]*F[0][0]*F[2][0] + a.d[5]*F[0][2]*F[2][0];

		cms /= pt.m_J;

		// LTE - Subtract away the stress term
		cms.d[0] -= s.xx();
		cms.d[9] -= s.yy();
		cms.d[20] -= s.zz();
		cms.d[6] -= s.xy();
		cms.d[18] -= s.yz();
		cms.d[15] -= s.xz();

		cms.d[27] -= s.xx();
		cms.d[2] -= s.yy();
		cms.d[14] -= s.zz();
		cms.d[22] -= s.xy();
		cms.d[11] -= s.yz();
		cms.d[25] -= s.xz();

		cms.d[44] -= s.xx();
		cms.d[35] -= s.yy();
		cms.d[5] -= s.zz();
		cms.d[43] -= s.xy();
		cms.d[30] -= s.yz();
		cms.d[38] -= s.xz();
		
		// LTE - Super-symmetrize to obtain the spatial elasticity tensor (this is where we obtain minor symmetry)
		c = cms.supersymm();
	}

	return c;
}
