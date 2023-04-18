/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



#include "stdafx.h"
#include "FERVEModel.h"
#include "FECore/FESolidDomain.h"
#include "FECore/FEElemElemList.h"
#include "FEBioMech/FEElasticMaterial.h"
#include "FEPeriodicBoundary1O.h"
#include "FECore/FEAnalysis.h"
#include "FECore/FELoadCurve.h"
#include "FEBioMech/FEBCPrescribedDeformation.h"
#include "FEBioMech/FESolidSolver2.h"
#include "FEBioMech/FEElasticSolidDomain.h"
#include <FECore/FEPeriodicLinearConstraint.h>
#include <FECore/FECube.h>
#include <FECore/FEPointFunction.h>
#include <FECore/FECoreKernel.h>

//-----------------------------------------------------------------------------
FERVEModel::FERVEModel()
{
	m_bctype = DISPLACEMENT;
}

//-----------------------------------------------------------------------------
FERVEModel::~FERVEModel()
{
}

//-----------------------------------------------------------------------------
void FERVEModel::SetParentModel(FEModel* fem)
{
	m_parentfem = fem;
}

//-----------------------------------------------------------------------------
// copy from the parent RVE
void FERVEModel::CopyFrom(FERVEModel& rve)
{
	// base class does most work
	FEModel::CopyFrom(rve);

	// copy the rest
	m_parentfem = rve.m_parentfem;
	m_bctype = rve.m_bctype;
	m_V0 = rve.m_V0;
	m_bb = rve.m_bb;
	m_BN = rve.m_BN;
}

//-----------------------------------------------------------------------------
bool FERVEModel::Init()
{
	// the RVE should only have one step
	if (Steps() != 1) return false;
	FEAnalysis* step = GetStep(0);

	// match the end time with the parent's
	FEAnalysis* parentStep = m_parentfem->GetStep(0);
	double tend = parentStep->m_ntime * parentStep->m_dt0;
	step->m_ntime = parentStep->m_ntime;
	step->m_dt0   = parentStep->m_dt0;
	step->m_tend  = tend;

	return FEModel::Init();
}

//-----------------------------------------------------------------------------
//! Initializes the RVE model and evaluates some useful quantities.
bool FERVEModel::InitRVE(int rveType, const char* szbc)
{
	// the RVE should only have one step
	if (Steps() != 1) return false;
	FEAnalysis* step = GetStep(0);

	// make sure the RVE problem doesn't output anything to a plot file
	step->SetPlotLevel(FE_PLOT_NEVER);

	// Center the RVE about the origin.
	// This also calculates the bounding box
	CenterRVE();

	// generate prescribed BCs
	// TODO: Make this part of the RVE definition
	m_bctype = rveType;
	if (m_bctype == DISPLACEMENT)
	{
		// find the boundary nodes
		if ((szbc) && (szbc[0] != 0))
		{
			// get the RVE mesh
			FEMesh& m = GetMesh();

			// find the node set that defines the corner nodes
			FENodeSet* pset = m.FindNodeSet(szbc);
			if (pset == 0) return false;

			// prep displacement BC's
			if (PrepDisplacementBC(pset) == false) return false;

			// tag all boundary nodes
			int NN = m.Nodes();
			m_BN.assign(NN, 0);
			FENodeSet& ns = *pset;
			for (int i=0; i<pset->Size(); ++i) m_BN[ns[i]] = 1; 
		}
		else 
		{
			// find all boundary nodes
			FindBoundaryNodes(m_BN);
			
			// create a (temporary) node set from the boundary nodes
			FEMesh& mesh = GetMesh();
			FENodeSet* set = new FENodeSet(this);

			int NN = mesh.Nodes();
			for (int i = 0; i<NN; ++i) if (m_BN[i] == 1) set->Add(i);

			// prep the displacement BCs
			if (PrepDisplacementBC(set) == false) return false;
		}
	}
	else if (m_bctype == PERIODIC_LC)
	{
		PrepPeriodicLC();
	}
/*	else if (m_bctype == PERIODIC_AL)
	{
		// prep periodic BC's
		if (PrepPeriodicBC(szbc) == false) return false;
	}
*/	else return false;

	// initialize base class
	if (FEModel::Init() == false) return false;

	// calculate intial RVE volume
	EvalInitialVolume();

	return true;
}

//-----------------------------------------------------------------------------
// This function sets up the periodic linear constraints boundary conditions.
// It assumes that the geometry is a cube.
bool FERVEModel::PrepPeriodicLC()
{
	// make sure there no BCs defined
	ClearBoundaryConditions();

	// user needs to define corner nodes
	// get the RVE mesh
	FEMesh& m = GetMesh();

	// Assuming it's a cube, build the surface, edge, and corner node data
	FECube cube;
	if (cube.Build(this) == false) return false;

	// tag all boundary nodes
	int NN = m.Nodes();
	const FENodeSet& bs = cube.GetBoundaryNodes();
	m_BN.resize(NN, 0);
	for (int i = 0; i<bs.Size(); ++i) m_BN[bs[i]] = 1;

	// now, build the linear constraints
	FEPeriodicLinearConstraint plc(this);
	plc.AddNodeSetPair(cube.GetSurface(0)->GetNodeList(), cube.GetSurface(1)->GetNodeList());
	plc.AddNodeSetPair(cube.GetSurface(2)->GetNodeList(), cube.GetSurface(3)->GetNodeList());
	plc.AddNodeSetPair(cube.GetSurface(4)->GetNodeList(), cube.GetSurface(5)->GetNodeList());
	plc.GenerateConstraints(this);

	// find the node set that defines the corner nodes
	FENodeSet* corners = const_cast<FENodeSet*>(&cube.GetCornerNodes());
	if (PrepDisplacementBC(corners) == false) return false;

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
	for (int n = 0; n < NN; ++n)
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
			int nf = el.Faces();
			for (int j=0; j<nf; ++j)
			{
				if (EEL.Neighbor(M, j) == 0)
				{
					// mark all nodes
					int nn = el.GetFace(j, fn);
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
bool FERVEModel::PrepDisplacementBC(FENodeSet* ns)
{
	// create a load curve
	FELoadCurve* plc = fecore_alloc(FELoadCurve, this);
	plc->Add(0.0, 0.0);
	plc->Add(1.0, 1.0);
	AddLoadController(plc);
	int NLC = LoadControllers() - 1;

	// clear all BCs
	ClearBoundaryConditions();

	// we create the prescribed deformation BC
	FEBCPrescribedDeformation* pdc = fecore_new<FEBCPrescribedDeformation>("prescribed deformation", this);
	AddBoundaryCondition(pdc);

	// assign the boundary nodes
	pdc->SetNodeSet(ns);

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
	int nc = SurfacePairConstraints();
	if (nc != 3) return false;
		
	for (int i=0; i<3; ++i)
	{
		FEPeriodicBoundary1O* pbc = dynamic_cast<FEPeriodicBoundary1O*>(SurfacePairConstraint(i));
		if (pbc == 0) return false;
	}

	// create a load curve
	FELoadCurve* plc = fecore_alloc(FELoadCurve, this);
	plc->Add(0.0, 0.0);
	plc->Add(1.0, 1.0);
	AddLoadController(plc);
	int NLC = LoadControllers() - 1;

	// create the DC's
	ClearBoundaryConditions();
	FEBCPrescribedDeformation* pdc = fecore_new<FEBCPrescribedDeformation>("prescribed deformation", this);
	AddBoundaryCondition(pdc);

	// assign nodes to BCs
	pdc->SetNodeSet(pset);

	// create the boundary node flags
	m_BN.assign(m.Nodes(), 0);
	int N = ns.Size();
	for (int i=0; i<N; ++i) m_BN[ns[i]] = 1;

	return true;
}

//-----------------------------------------------------------------------------
void FERVEModel::Update(const mat3d& F)
{
	// get the mesh
	FEMesh& m = GetMesh();

	// assign new DC's for the boundary nodes
	FEBCPrescribedDeformation& dc = dynamic_cast<FEBCPrescribedDeformation&>(*BoundaryCondition(0));
	dc.SetDeformationGradient(F);

	if (m_bctype == FERVEModel::PERIODIC_AL)
	{
		// loop over periodic boundaries
		for (int i = 0; i<3; ++i)
		{
			FEPeriodicBoundary1O* pc = dynamic_cast<FEPeriodicBoundary1O*>(SurfacePairConstraint(i));
			assert(pc);

			pc->m_Fmacro = F;
		}
	}
}

//-----------------------------------------------------------------------------
void FERVEModel::ScaleGeometry(double scale)
{
	// get the mesh
	FEMesh& mesh = GetMesh();

	// calculate the center of mass first
	vec3d rc(0,0,0);
	for (int i = 0; i<mesh.Nodes(); ++i)
	{
		rc += mesh.Node(i).m_r0;
	}
	rc /= (double) mesh.Nodes();

	// scale the nodal positions around the center
	for (int i = 0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		vec3d r0 = node.m_r0;
		node.m_r0 = rc + (r0 - rc)*scale;

		vec3d rt = node.m_rt;
		node.m_rt = rc + (rt - rc)*scale;
	}
}

//-----------------------------------------------------------------------------
//! return current volume (calculated each time)
double FERVEModel::CurrentVolume()
{
	// get the RVE mesh
	FEMesh& m = GetMesh();

	double V = 0.0;
	for (int i = 0; i<m.Domains(); ++i)
	{
		FESolidDomain& dom = dynamic_cast<FESolidDomain&>(m.Domain(i));
		int NE = dom.Elements();
		for (int j = 0; j<NE; ++j)
		{
			FESolidElement& el = dom.Element(j);

			int ni = el.GaussPoints();
			int ne = el.Nodes();
			double* w = el.GaussWeights();
			for (int n = 0; n<ni; ++n)
			{
				FEMaterialPoint& pt = *el.GetMaterialPoint(n);
				FEElasticMaterialPoint& ep = *pt.ExtractData<FEElasticMaterialPoint>();

				// calculate Jacobian
				double Jn = dom.detJt(el, n);

				// add it all up
				V += w[n] * Jn;
			}
		}
	}

	return V;
}

//-----------------------------------------------------------------------------
mat3ds FERVEModel::StressAverage(mat3d& F, FEMaterialPoint& mp)
{
	// rewind the RCI
	RCI_Rewind();

	// update the BC's
	Update(F);

	// match the time step
	FEModel& parent = *m_parentfem;
	FETimeInfo& ti = parent.GetTime();
	SetCurrentTimeStep(ti.timeIncrement);

	// advance the RVE solution
	bool bret = RCI_Advance();

	assert(ti.currentTime == GetCurrentTime());

	// make sure it converged
	if (bret == false) throw FEMultiScaleException(-1, -1);

	// calculate and return the (Cuachy) stress average
	return StressAverage(mp);
}

//-----------------------------------------------------------------------------
//! Calculate the stress average
mat3ds FERVEModel::StressAverage(FEMaterialPoint& mp)
{
	// get the RVE mesh
	FEMesh& m = GetMesh();

	mat3ds T; T.zero();

	for (int i=0; i<m.Domains(); ++i)
	{
		FESolidDomain& dom = dynamic_cast<FESolidDomain&>(m.Domain(i));
		int NE = dom.Elements();
		for (int j=0; j<NE; ++j)
		{
			FESolidElement& el = dom.Element(j);

			int ni = el.GaussPoints();
			int ne = el.Nodes();
			double* w = el.GaussWeights();
			for (int n=0; n<ni; ++n)
			{
				FEMaterialPoint& pt = *el.GetMaterialPoint(n);
				FEElasticMaterialPoint& ep = *pt.ExtractData<FEElasticMaterialPoint>();

				// calculate Jacobian
				double Jn = dom.detJt(el, n);

				// add it all up
				T += ep.m_s*(w[n]*Jn);
			}
		}
	}

/*	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	double J = pt.m_J;
	double V0 = InitialVolume();
	return T / (J*V0);
*/
	return T / CurrentVolume();
}

/*
//-----------------------------------------------------------------------------
//! Calculate the stress average
mat3ds FERVEModel::StressAverage(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	double J = pt.m_J;

	// get the RVE mesh
	FEMesh& m = GetMesh();

	mat3d T; T.zero();

	// for periodic BC's we take the reaction forces directly from the periodic constraints
	// TODO: Figure out a way to store all the reaction forces on the nodes
	//       That way, we don't need to do this anymore.
	if (m_bctype == FERVEModel::PERIODIC_AL)
	{
		// get the reaction for from the periodic constraints
		for (int i = 0; i<3; ++i)
		{
			FEPeriodicBoundary1O* pbc = dynamic_cast<FEPeriodicBoundary1O*>(SurfacePairInteraction(i));
			assert(pbc);
			FEPeriodicSurface& ss = pbc->m_ss;
			int N = ss.Nodes();
			for (int i = 0; i<N; ++i)
			{
				FENode& node = ss.Node(i);
				vec3d f = ss.m_Fr[i];

				// We multiply by two since the reaction forces are only stored at the primary surface 
				// and we also need to sum over the secondary nodes (NOTE: should I figure out a way to 
				// store the reaction forces on the secondary nodes as well?)
				T += (f & node.m_rt)*2.0;
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

	if (m_bctype != FERVEModel::PERIODIC_LC)
	{
		// calculate the averaged stress
		// TODO: This might be more efficient if we keep a list of boundary nodes
		int NN = m.Nodes();
		for (int j = 0; j<NN; ++j)
		{
			if (IsBoundaryNode(j))
			{
				FENode& n = m.Node(j);
				vec3d f;
				f.x = R[-n.m_ID[dof_X] - 2];
				f.y = R[-n.m_ID[dof_Y] - 2];
				f.z = R[-n.m_ID[dof_Z] - 2];
				T += f & n.m_rt;
			}
		}
	}
	else
	{
		for (int i = 0; i<m_FN.size(); ++i)
		{
			FENode& n = m.Node(m_FN[i]);
			vec3d f;
			f.x = R[-n.m_ID[dof_X] - 2];
			f.y = R[-n.m_ID[dof_Y] - 2];
			f.z = R[-n.m_ID[dof_Z] - 2];
			T += f & n.m_rt;
		}
	}

	double V0 = InitialVolume();
	return T.sym() / (J*V0);
}
*/


//-----------------------------------------------------------------------------
//! Calculate the average stiffness from the RVE solution.
tens4ds FERVEModel::StiffnessAverage(FEMaterialPoint &mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// get the mesh
	FEMesh& m = GetMesh();

	// the element's stiffness matrix
	matrix ke;

	// element's residual
	vector<double> fe;

	// calculate the center point
	vec3d rc(0, 0, 0);
	for (int k = 0; k<m.Nodes(); ++k) rc += m.Node(k).m_rt;
	rc /= (double)m.Nodes();

	// elasticity tensor
	tens4ds c(0.);

	// calculate the stiffness matrix and residual
	for (int k = 0; k<m.Domains(); ++k)
	{
		FEElasticSolidDomain& bd = dynamic_cast<FEElasticSolidDomain&>(m.Domain(k));
		int NS = bd.Elements();
		for (int n = 0; n<NS; ++n)
		{
			FESolidElement& e = bd.Element(n);

			// create the element's stiffness matrix
			int ne = e.Nodes();
			int ndof = 3 * ne;
			ke.resize(ndof, ndof);
			ke.zero();

			// calculate the element's stiffness matrix
			bd.ElementStiffness(GetTime(), n, ke);

			// create the element's residual
			fe.assign(ndof, 0);

			// calculate the element's residual
			bd.ElementInternalForce(e, fe);

			// loop over the element's nodes
			for (int i = 0; i<ne; ++i)
			{
				FENode& ni = m.Node(e.m_node[i]);
				for (int j = 0; j<ne; ++j)
				{
					FENode& nj = m.Node(e.m_node[j]);
					if (IsBoundaryNode(e.m_node[i]) && IsBoundaryNode(e.m_node[j]))
					{
						// both nodes are boundary nodes
						// so grab the element's submatrix
						mat3d K;
						ke.get(3 * i, 3 * j, K);

						// get the nodal positions
						vec3d ri = ni.m_rt - rc;
						vec3d rj = nj.m_rt - rc;

						// create the elasticity tensor
						c += dyad4s(ri, K, rj);
					}
				}
			}
		}
	}

	// divide by volume
/*	double V0 = InitialVolume();
	double Vi = 1.0 / (pt.m_J * V0);
*/
	double Vi = 1.0 / CurrentVolume();
	c *= Vi;

	return c;
}

// this function is hidden
bool FERVEModel::Solve() { assert(false); return false; }
