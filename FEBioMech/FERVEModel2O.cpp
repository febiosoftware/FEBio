#include "stdafx.h"
#include "FERVEModel2O.h"
#include "FECore/FESolidDomain.h"
#include "FECore/FEElemElemList.h"
#include "FECore/BC.h"
#include "FEElasticMaterial.h"
#include "FEPeriodicBoundary2O.h"
#include "FECore/FEAnalysis.h"
#include "FECore/LoadCurve.h"
#include "FESolidSolver2.h"
#include "FEElasticSolidDomain.h"

//-----------------------------------------------------------------------------
FERVEModel2O::FERVEModel2O()
{
	m_bperiodic = false;

	// set the pardiso solver as default
	m_nsolver = PARDISO_SOLVER;
}

//-----------------------------------------------------------------------------
FERVEModel2O::~FERVEModel2O()
{
}

//-----------------------------------------------------------------------------
//! Initializes the RVE model and evaluates some useful quantities.
bool FERVEModel2O::InitRVE(bool bperiodic, const char* szbc)
{
	// make sure the RVE problem doesn't output anything to a plot file
	GetCurrentStep()->SetPlotLevel(FE_PLOT_NEVER);

	// Center the RVE about the origin.
	// This also calculates the bounding box
	CenterRVE();

	// generate prescribed BCs
	// TODO: Make this part of the RVE definition
	m_bperiodic = bperiodic;
	if (bperiodic == false)
	{
		// find the boundary nodes
		FindBoundaryNodes(m_BN);

		// prep displacement BC's
		if (PrepDisplacementBC() == false) return false;
	}
	else
	{
		// prep periodic BC's
		if (PrepPeriodicBC(szbc) == false) return false;
	}

	// initialize base class
	if (FEModel::Init() == false) return false;

	// calculate intial RVE volume
	EvalInitialVolume();

	return true;
}

//-----------------------------------------------------------------------------
//! Evaluates the initial volume of the RVE model.
//! This is called from FERVEModel2O::Init.
void FERVEModel2O::EvalInitialVolume()
{
	m_V0 = 0;
	FEMesh& m = GetMesh();
	for (int k=0; k<m.Domains(); ++k)
	{
		FESolidDomain& dom = static_cast<FESolidDomain&>(m.Domain(k));
		for (int i=0; i<dom.Elements(); ++i)
		{
			FESolidElement& el = dom.Element(i);
			int nint = el.GaussPoints();
			double* w = el.GaussWeights();
			double ve = 0;
			for (int n=0; n<nint; ++n)
			{
				FEElasticMaterialPoint& pt = *el.GetMaterialPoint(n)->ExtractData<FEElasticMaterialPoint>();
				double J = dom.detJ0(el, n);
				ve += J*w[n];
			}
			m_V0 += ve;
		}
	}
}

//-----------------------------------------------------------------------------
//! Centers the RVE around the origin.
void FERVEModel2O::CenterRVE()
{
	FEMesh& mesh = GetMesh();
	FENode& node = mesh.Node(0);

	// setup bounding box
	FEBoundingBox box(node.m_r0);
	const int NN = mesh.Nodes();
	for (int i=1; i<NN; ++i)
	{
		FENode& node = mesh.Node(i);
		box.add(node.m_r0);
	}

	// get the center
	vec3d c = box.center();
	
	// recenter the RVE about the origin
	for (int n=0; n<NN; ++n)
	{
		FENode& node = mesh.Node(n);
		node.m_r0 -= c;
		node.m_rt = node.m_r0;
	}

	// adjust bounding box
	m_bb.translate(-c);
}

//-----------------------------------------------------------------------------
//! Find the boundary nodes of the RVE model
void FERVEModel2O::FindBoundaryNodes(vector<int>& BN)
{
	// first we need to find all the boundary nodes
	FEMesh& m = GetMesh();
	int NN = m.Nodes();
	BN.assign(NN, 0);

	// create the element-element list
	FEElemElemList EEL;
	EEL.Create(&m);

	double wx = m_bb.width()*0.5;
	double wy = m_bb.height()*0.5;
	double wz = m_bb.depth()*0.5;

	// use the E-E list to tag all exterior nodes
	int fn[FEElement::MAX_NODES], M = 0;
	for (int k=0; k<m.Domains(); ++k)
	{
		FEDomain& dom = m.Domain(k);
		for (int i=0; i<dom.Elements(); ++i, ++M)
		{
			FEElement& el = dom.ElementRef(i);
			int nf = m.Faces(el);
			for (int j=0; j<nf; ++j)
			{
				if (EEL.Neighbor(M, j) == 0)
				{
					// mark all nodes
					int nn = m.GetFace(el, j, fn);
					for (int k=0; k<nn; ++k)
					{
						FENode& node = m.Node(fn[k]);
						
						if (fabs(node.m_r0.x) >= 0.999*wx) BN[fn[k]] = 1;
						if (fabs(node.m_r0.y) >= 0.999*wy) BN[fn[k]] = 1;
						if (fabs(node.m_r0.z) >= 0.999*wz) BN[fn[k]] = 1;
					}
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------
// Setup the displacement boundary conditions.
bool FERVEModel2O::PrepDisplacementBC()
{
	FEMesh& m = GetMesh();
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
	AddLoadCurve(plc);
	int NLC = LoadCurves() - 1;

	// clear all BCs
	ClearBCs();

	// we create three DCs, one for each displacement dof
	FEPrescribedBC* pdc[3] = {0};
	pdc[0] = new FEPrescribedBC(this); pdc[0]->SetDOF(0).SetScale(1.0, NLC); AddPrescribedBC(pdc[0]);
	pdc[1] = new FEPrescribedBC(this); pdc[1]->SetDOF(1).SetScale(1.0, NLC); AddPrescribedBC(pdc[1]);
	pdc[2] = new FEPrescribedBC(this); pdc[2]->SetDOF(2).SetScale(1.0, NLC); AddPrescribedBC(pdc[2]);

	// assign the boundary nodes
	for (i=0; i<N; ++i)
		if (m_BN[i] == 1)
		{
			pdc[0]->AddNode(i, 0.0);
			pdc[1]->AddNode(i, 0.0);
			pdc[2]->AddNode(i, 0.0);
		}

	return true;
}

//-----------------------------------------------------------------------------
bool FERVEModel2O::PrepPeriodicBC(const char* szbc)
{
	// make sure the node set is valid
	if ((szbc==0)||(szbc[0]==0)) return false;

	// get the RVE mesh
	FEMesh& m = GetMesh();

	// find the node set that defines the corner nodes
	FENodeSet* pset = m.FindNodeSet(szbc);
	if (pset == 0) return false;
	FENodeSet& ns = *pset;

	// check the periodic constraints
	int nc = SurfacePairInteractions();
	if (nc != 3) return false;
		
	for (int i=0; i<3; ++i)
	{
		FEPeriodicBoundary2O* pbc = dynamic_cast<FEPeriodicBoundary2O*>(SurfacePairInteraction(i));
		if (pbc == 0) return false;
	}

	// create a load curve
	FELoadCurve* plc = new FELoadCurve;
	plc->SetInterpolation(FELoadCurve::LINEAR);
	plc->Add(0.0, 0.0);
	plc->Add(1.0, 1.0);
	AddLoadCurve(plc);
	int NLC = LoadCurves() - 1;

	// create the DC's
	ClearBCs();
	FEPrescribedBC* pdc[3] = {0};
	pdc[0] = new FEPrescribedBC(this); pdc[0]->SetDOF(0).SetScale(1.0, NLC); AddPrescribedBC(pdc[0]);
	pdc[1] = new FEPrescribedBC(this); pdc[1]->SetDOF(1).SetScale(1.0, NLC); AddPrescribedBC(pdc[1]);
	pdc[2] = new FEPrescribedBC(this); pdc[2]->SetDOF(2).SetScale(1.0, NLC); AddPrescribedBC(pdc[2]);

	// assign nodes to BCs
	pdc[0]->AddNodes(ns, 0.0);
	pdc[1]->AddNodes(ns, 0.0);
	pdc[2]->AddNodes(ns, 0.0);

	// create the boundary node flags
	m_BN.assign(m.Nodes(), 0);
	int N = ns.size();
	for (int i=0; i<N; ++i) m_BN[ns[i]] = 1;

	return true;
}

//=============================================================================
FEMicroModel2O::FEMicroModel2O()
{
	m_bperiodic = false;
	m_V0 = 0.0;
	m_BN = 0;
}

//-----------------------------------------------------------------------------
FEMicroModel2O::~FEMicroModel2O()
{
}

//-----------------------------------------------------------------------------
bool FEMicroModel2O::Init(FERVEModel2O& rve)
{
	// copy the master RVE
	CopyFrom(rve);

	// initialize model
	if (FEModel::Init() == false) return false;

	// copy some stuff from the master RVE model
	m_V0 = rve.InitialVolume();
	m_bperiodic = rve.IsPeriodic();
	m_BN = &rve.BoundaryList();
	if (m_BN == 0) return false;

	return true;
}

//-----------------------------------------------------------------------------
//! Solve the RVE model
bool FEMicroModel2O::Solve(const mat3d& F, const tens3drs& G)
{
	// update boundary conditions
	UpdateBC(F, G);

	// solve the model
	return FEModel::Solve();
}

//-----------------------------------------------------------------------------
//! Assign the prescribed displacement to the boundary nodes.
void FEMicroModel2O::UpdateBC(const mat3d& F, const tens3drs& G)
{
	// get the mesh
	FEMesh& m = GetMesh();

	// assign new DC's for the boundary nodes
	FEPrescribedBC& dx = *PrescribedBC(0);
	FEPrescribedBC& dy = *PrescribedBC(1);
	FEPrescribedBC& dz = *PrescribedBC(2);

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
			FEPeriodicBoundary2O* pc = dynamic_cast<FEPeriodicBoundary2O*>(SurfacePairInteraction(i));
			assert(pc);
			pc->m_Fmacro = F;
			pc->m_Gmacro = G;
		}
	}
}

//-----------------------------------------------------------------------------
//! Calculate the average stress from the RVE solution.
mat3d FEMicroModel2O::AveragedStressPK1(FEMaterialPoint &mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	mat3d F = pt.m_F;
	double J = pt.m_J;
	
	// get the RVE mesh
	FEMesh& m = GetMesh();

	mat3d PK1; PK1.zero();

	// for periodic BC's we take the reaction forces directly from the periodic constraints
	if (m_bperiodic)
	{
		// get the reaction for from the periodic constraints
		for (int i=0; i<3; ++i)
		{
			FEPeriodicBoundary2O* pbc = dynamic_cast<FEPeriodicBoundary2O*>(SurfacePairInteraction(i));
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
	const int dof_X = GetDOFIndex("x");
	const int dof_Y = GetDOFIndex("y");
	const int dof_Z = GetDOFIndex("z");
	FEAnalysis* pstep = GetCurrentStep();
	FESolidSolver2* ps = dynamic_cast<FESolidSolver2*>(pstep->GetFESolver());
	assert(ps);
	vector<double>& R = ps->m_Fr;
	FEPrescribedBC& dc = *PrescribedBC(0);
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

	return PK1 / m_V0;
}

//-----------------------------------------------------------------------------
//! Calculate the average stress from the RVE solution.
void FEMicroModel2O::AveragedStress2O(mat3d& Pa, tens3drs& Qa)
{
	// get the RVE mesh
	FEMesh& m = GetMesh();

	Pa.zero();
	Qa.zero();

	// for periodic BC's we take the reaction forces directly from the periodic constraints
	if (m_bperiodic)
	{
		// get the reaction for from the periodic constraints
		for (int i=0; i<3; ++i)
		{
			FEPeriodicBoundary2O* pbc = dynamic_cast<FEPeriodicBoundary2O*>(SurfacePairInteraction(i));
			assert(pbc);
			FEPeriodicSurface& ss = pbc->m_ss;
			int N = ss.Nodes();
			for (int i=0; i<N; ++i)
			{
				FENode& node = ss.Node(i);
				vec3d f = ss.m_Fr[i];
				const vec3d& X = node.m_r0;
				
				// We multiply by two since the reaction forces are only stored at the slave surface 
				// and we also need to sum over the master nodes (NOTE: should I figure out a way to 
				// store the reaction forces on the master nodes as well?)
				Pa += (f & X)*2.0;

				Qa += dyad3rs(f, X)*2.0;
			}
		}
	}

	// get the reaction force vector from the solid solver
	// (We also need to do this for the periodic BC, since at the prescribed nodes,
	// the contact forces will be zero).
	FEPrescribedBC& dx = *PrescribedBC(0);
	FEPrescribedBC& dy = *PrescribedBC(1);
	FEPrescribedBC& dz = *PrescribedBC(2);

	const int dof_X = GetDOFIndex("x");
	const int dof_Y = GetDOFIndex("y");
	const int dof_Z = GetDOFIndex("z");
	FEAnalysis* pstep = GetCurrentStep();
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

		const vec3d& X = n.m_r0;

		Pa += (f & X);

		Qa += dyad3rs(f, X);
	}

	// divide by volume
	Pa /= m_V0;
	Qa /= 2*m_V0;
}
/*
//-----------------------------------------------------------------------------
//! Calculate the average stress from the RVE solution.
void FEMicroModel2O::AveragedStress2OPK1(FEMaterialPoint &mp, mat3d &PK1a, tens3drs &QK1a)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	mat3d F = pt.m_F;
	double J = pt.m_J;

	// get the RVE mesh
	FEMesh& m = GetMesh();

	mat3d PK1; PK1.zero();
	tens3drs QK1; QK1.zero();

	// for periodic BC's we take the reaction forces directly from the periodic constraints
	if (m_bperiodic)
	{
		// get the reaction for from the periodic constraints
		for (int i=0; i<3; ++i)
		{
			FEPeriodicBoundary2O* pbc = dynamic_cast<FEPeriodicBoundary2O*>(SurfacePairInteraction(i));
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
	FEPrescribedBC& dx = *PrescribedBC(0);
	FEPrescribedBC& dy = *PrescribedBC(1);
	FEPrescribedBC& dz = *PrescribedBC(2);

	const int dof_X = GetDOFIndex("x");
	const int dof_Y = GetDOFIndex("y");
	const int dof_Z = GetDOFIndex("z");
	FEAnalysis* pstep = GetCurrentStep();
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

	PK1a = PK1 / m_V0;
	QK1a = QK1 / (2*m_V0);
}

//-----------------------------------------------------------------------------
//! Calculate the average stress from the RVE solution.
void FEMicroModel2O::AveragedStress2OPK2(FEMaterialPoint &mp, mat3ds &Sa, tens3ds &Ta)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	mat3d F = pt.m_F;
	double J = pt.m_J;
	mat3d Finv = F.inverse();

	// get the RVE mesh
	FEMesh& m = GetMesh();

	mat3d S; S.zero();
	tens3ds T; T.zero();

	// for periodic BC's we take the reaction forces directly from the periodic constraints
	if (m_bperiodic)
	{
		// get the reaction for from the periodic constraints
		for (int i=0; i<3; ++i)
		{
			FEPeriodicBoundary2O* pbc = dynamic_cast<FEPeriodicBoundary2O*>(SurfacePairInteraction(i));
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
	FEPrescribedBC& dx = *PrescribedBC(0);
	FEPrescribedBC& dy = *PrescribedBC(1);
	FEPrescribedBC& dz = *PrescribedBC(2);

	const int dof_X = GetDOFIndex("x");
	const int dof_Y = GetDOFIndex("y");
	const int dof_Z = GetDOFIndex("z");
	FEAnalysis* pstep = GetCurrentStep();
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

	Sa = S.sym() / m_V0;
	Ta = T / (2*m_V0);
}
*/

//-----------------------------------------------------------------------------
//! Calculate the average stiffness from the RVE solution.
void FEMicroModel2O::AveragedStiffness(FEMaterialPoint &mp, tens4d& C, tens5d& L, tens5d& H, tens6d& J)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// get the mesh
	FEMesh& m = GetMesh();

	// get the solver
	FEAnalysis* pstep = GetCurrentStep();
	FESolidSolver2* ps = dynamic_cast<FESolidSolver2*>(pstep->GetFESolver());

	// the element's stiffness matrix
	matrix ke;

	// element's residual
	vector<double> fe;

	// calculate the center point
	vec3d rc(0,0,0);
	for (int k=0; k<m.Nodes(); ++k) rc += m.Node(k).m_r0;
	rc /= (double) m.Nodes();

	// zero the stiffness components
	C.zero();
	L.zero();
	H.zero();
	J.zero();

	// calculate the stiffness matrix and residual
	for (int nd=0; nd<m.Domains(); ++nd)
	{
		FEElasticSolidDomain& bd = dynamic_cast<FEElasticSolidDomain&>(m.Domain(nd));
		int NEL = bd.Elements();
		for (int ne=0; ne<NEL; ++ne)
		{
			FESolidElement& el = bd.Element(ne);

			// create the element's stiffness matrix
			int neln = el.Nodes();
			int ndof = 3*neln;
			ke.resize(ndof, ndof);
			ke.zero();

			// calculate the element's stiffness matrix
			bd.ElementStiffness(GetTime(), ne, ke);

			// create the element's residual
			fe.assign(ndof, 0);

			// calculate the element's residual
			bd.ElementInternalForce(el, fe);

			// loop over the element's nodes
			for (int a=0; a<neln; ++a)
			{
				FENode& na = m.Node(el.m_node[a]);
				for (int b=0; b<neln; ++b)
				{
					FENode& nb = m.Node(el.m_node[b]);
					if (IsBoundaryNode(el.m_node[a]) && IsBoundaryNode(el.m_node[b]))
					{
						// both nodes are boundary nodes
						// so grab the element's submatrix
						mat3d K;
						ke.get(3*a, 3*b, K);

						// get the nodal positions relative to the center
						vec3d ra = na.m_r0 - rc;
						vec3d rb = nb.m_r0 - rc;
						
						double Ra[3] = { ra.x, ra.y, ra.z };
						double Rb[3] = { rb.x, rb.y, rb.z };

						// create the elasticity tensors
						for (int i=0; i<3; ++i)
							for (int j=0; j<3; ++j)
								for (int k=0; k<3; ++k)
									for (int l=0; l<3; ++l)
									{
										C(i, j, k, l) += K[i][k]*Rb[l]*Ra[j];
									}

						for (int i=0; i<3; ++i)
							for (int j=0; j<3; ++j)
								for (int k=0; k<3; ++k)
									for (int l=0; l<3; ++l)
										for (int m=0; m<3; ++m)
										{
											L(i, j, k, l, m) += K[i][k]*Rb[l]*Rb[m]*Ra[j];
										}

						for (int i=0; i<3; ++i)
							for (int j=0; j<3; ++j)
								for (int k=0; k<3; ++k)
									for (int l=0; l<3; ++l)
										for (int m=0; m<3; ++m)
										{
											H(i, j, k, l, m) += K[i][l]*Rb[m]*Ra[j]*Ra[k];
										}

						for (int i=0; i<3; ++i)
							for (int j=0; j<3; ++j)
								for (int k=0; k<3; ++k)
									for (int l=0; l<3; ++l)
										for (int m=0; m<3; ++m)
											for (int n=0; n<3; ++n)
											{
												J(i, j, k, l, m, n) += K[i][l]*Rb[m]*Rb[n]*Ra[j]*Ra[k];
											}

					}
				}
			}
		}
	}

	// divide by volume
	C /= m_V0;
	L /= m_V0;
	H /= 2.0*m_V0;
	J /= 2.0*m_V0;
}
