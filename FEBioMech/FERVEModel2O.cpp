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
	FEBoundingBox box(node.m_r0, node.m_r0);
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
void FEMicroModel2O::AveragedStress2O(FEMaterialPoint &mp, mat3ds &sa, tens3ds &taua)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	mat3d F = pt.m_F;
	double J = pt.m_J;

	// get the RVE mesh
	FEMesh& m = GetMesh();

	mat3d s; s.zero();
	tens3ds tau; tau.zero();

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
				s += (f & node.m_rt)*2.0;

				vec3d x; x = node.m_rt;
		
				tau += dyad3s(x, f, x)*2.0;
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
				
		s += (f & n.m_rt);

		vec3d x = n.m_rt;
		
		tau += dyad3s(x, f, x);
	}

	sa = s.sym() / (J*m_V0);
	taua = tau / (2*J*m_V0);
}

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

//-----------------------------------------------------------------------------
//! Calculate the average stiffness from the RVE solution.
void FEMicroModel2O::AveragedStiffness(FEMaterialPoint &mp, tens4ds& c, tens5ds& d, tens6ds& e)
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
			bd.ElementStiffness(GetTime(), n, ke);

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
					if (IsBoundaryNode(el.m_node[i]) && IsBoundaryNode(el.m_node[j]))
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
	c = tens4ds(D)/(pt.m_J * m_V0);
	d = d/(2.*pt.m_J * m_V0);
	e = e/(4.*pt.m_J * m_V0);
}
