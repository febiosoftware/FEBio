#include "stdafx.h"
#include "FERVEModel.h"
#include "FECore/FESolidDomain.h"
#include "FECore/FEElemElemList.h"
#include "FEElasticMaterial.h"
#include "FEPeriodicBoundary1O.h"

//-----------------------------------------------------------------------------
FERVEModel::FERVEModel()
{
	m_bperiodic = false;

	// set the pardiso solver as default
	m_nsolver = PARDISO_SOLVER;
}

//-----------------------------------------------------------------------------
FERVEModel::~FERVEModel()
{
}

//-----------------------------------------------------------------------------
//! Initializes the RVE model and evaluates some useful quantities.
bool FERVEModel::InitRVE(bool bperiodic, const char* szbc)
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
//! This is called from FERVEModel::Init.
void FERVEModel::EvalInitialVolume()
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
				double J = dom.detJt(el, n);
				ve += J*w[n];
			}
			m_V0 += ve;
		}
	}
}

//-----------------------------------------------------------------------------
//! Centers the RVE around the origin.
void FERVEModel::CenterRVE()
{
	FEMesh& mesh = GetMesh();
	FENode& node = mesh.Node(0);

	// setup bounding box
	FE_BOUNDING_BOX box;
	box.r0 = node.m_r0;
	box.r1 = box.r0;
	const int NN = mesh.Nodes();
	for (int i=1; i<NN; ++i)
	{
		FENode& node = mesh.Node(i);

		const vec3d& ri = node.m_r0;
		box += ri;
	}

	// get the center
	vec3d c = box.center();
	
	// recenter the RVE about the origin
	for (int n = 0; n < NN; ++n)
	{
		FENode& node = mesh.Node(n);
		node.m_r0 -= c;
		node.m_rt = node.m_r0;
	}

	// adjust bounding box
	m_bb.r0 = box.r0 - c;
	m_bb.r1 = box.r1 - c;
}

//-----------------------------------------------------------------------------
//! Find the boundary nodes of the RVE model
void FERVEModel::FindBoundaryNodes(vector<int>& BN)
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
bool FERVEModel::PrepDisplacementBC()
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
	pdc[0] = new FEPrescribedBC(this); pdc[0]->SetDOF(0).SetLoadCurveIndex(NLC).SetScale(1.0); AddPrescribedBC(pdc[0]);
	pdc[1] = new FEPrescribedBC(this); pdc[1]->SetDOF(1).SetLoadCurveIndex(NLC).SetScale(1.0); AddPrescribedBC(pdc[1]);
	pdc[2] = new FEPrescribedBC(this); pdc[2]->SetDOF(2).SetLoadCurveIndex(NLC).SetScale(1.0); AddPrescribedBC(pdc[2]);

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
bool FERVEModel::PrepPeriodicBC(const char* szbc)
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
		FEPeriodicBoundary1O* pbc = dynamic_cast<FEPeriodicBoundary1O*>(SurfacePairInteraction(i));
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
	pdc[0] = new FEPrescribedBC(this); pdc[0]->SetDOF(0).SetLoadCurveIndex(NLC).SetScale(1.0); AddPrescribedBC(pdc[0]);
	pdc[1] = new FEPrescribedBC(this); pdc[1]->SetDOF(1).SetLoadCurveIndex(NLC).SetScale(1.0); AddPrescribedBC(pdc[1]);
	pdc[2] = new FEPrescribedBC(this); pdc[2]->SetDOF(2).SetLoadCurveIndex(NLC).SetScale(1.0); AddPrescribedBC(pdc[2]);

	m_BN.assign(m.Nodes(), 0);

	int N = ns.size();
	for (int i=0; i<N; ++i)
	{
		m_BN[ns[i]] = 1;
		pdc[0]->AddNode(ns[i], 0.0);
		pdc[1]->AddNode(ns[i], 0.0);
		pdc[2]->AddNode(ns[i], 0.0);
	}

	return true;
}
