#include "stdafx.h"
#include "FEMicroMaterial2O.h"
#include "FECore/FEElemElemList.h"
#include "FECore/log.h"
#include "FESolidSolver.h"
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
	m_inf_str.zero();
	m_inf_str_grad.zero();
	m_E.zero();
	m_H.zero();
	m_e.zero();
	m_h.zero();
	m_energy_diff = 0.;

	m_Ca.zero();
	m_Da.zero();
	m_Ea.zero();
}

//-----------------------------------------------------------------------------
//! Initialize material point data
void FEMicroMaterialPoint2O::Init(bool bflag)
{
}

//-----------------------------------------------------------------------------
//! create a shallow copy
FEMaterialPoint* FEMicroMaterialPoint2O::Copy()
{
	FEMicroMaterialPoint2O* pt = new FEMicroMaterialPoint2O(m_pt?m_pt->Copy():0);
	pt->m_tau = m_tau;
	pt->m_G = m_G;
	pt->m_Ca = m_Ca;
	pt->m_Da = m_Da;
	pt->m_Ea = m_Ea;
	return pt;
}

//-----------------------------------------------------------------------------
//! serialize material point data
void FEMicroMaterialPoint2O::Serialize(DumpFile& ar)
{
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

//-----------------------------------------------------------------------------
//! stream material point data
void FEMicroMaterialPoint2O::ShallowCopy(DumpStream& dmp, bool bsave)
{
	if (bsave)
	{
		dmp << m_Ca;
		dmp << m_Da;
		dmp << m_Ea;
	}
	else
	{
		dmp >> m_Ca;
		dmp >> m_Da;
		dmp >> m_Ea;
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
FEMicroMaterial2O::FEMicroMaterial2O(FEModel* pfem) : FEElasticMaterial(pfem)
{
	m_brve = false;

	// initialize parameters
	m_szrve[0] = 0;
	m_szbc[0] = 0;
	m_bperiodic = false;

	m_bb_x = 0.; m_bb_y = 0.; m_bb_z = 0.;
	m_num_ext_node = 0;
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
void FEMicroMaterial2O::Init()
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
bool FEMicroMaterial2O::PrepRVE()
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
void FEMicroMaterial2O::FindBoundaryNodes()
{
	// first we need to find all the boundary nodes
	FEMesh& m = m_rve.GetMesh();
	int N = m.Nodes();
	m_BN.assign(N, 0);

	// create the element-element list
	FEElemElemList EEL;
	EEL.Create(&m);

	// LTE - Find the initial bounding box and center of the RVE
	double xmin = 0.; double ymin = 0.; double zmin = 0.;
	double xmax = 0.; double ymax = 0.; double zmax = 0.;
	double centx = 0.; double centy = 0.; double centz = 0.;

	FENode& node = m.Node(0);
	xmin = node.m_r0.x; xmax = node.m_r0.x; 
	ymin = node.m_r0.y; ymax = node.m_r0.y; 
	zmin = node.m_r0.z; zmax = node.m_r0.z; 

	for (int n = 1; n < N; ++n){
		FENode& node = m.Node(n);

		if (node.m_r0.x >= xmax) xmax = node.m_r0.x;
		if (node.m_r0.x <= xmin) xmin = node.m_r0.x;
		if (node.m_r0.y >= ymax) ymax = node.m_r0.y;
		if (node.m_r0.y <= ymin) ymin = node.m_r0.y;
		if (node.m_r0.z >= zmax) zmax = node.m_r0.z;
		if (node.m_r0.z <= zmin) zmin = node.m_r0.z;
	}

	centx = xmin + (xmax - xmin)/2; centy = ymin + (ymax - ymin)/2; centz = zmin + (zmax - zmin)/2;
	
	// LTE - Recenter the RVE about the origin
	for (int n = 0; n < N; ++n){
		FENode& node = m.Node(n);
		node.m_r0.x -= centx; node.m_r0.y -= centy; node.m_r0.z -= centz;
		node.m_rt.x = node.m_r0.x; node.m_rt.y = node.m_r0.y; node.m_rt.z = node.m_r0.z;
	}
	
	// LTE - Find the bounding box for the RVE (+/- Lx, Ly, Lz)
	m_bb_x = xmax - centx; m_bb_y = ymax - centy; m_bb_z = zmax - centz;
	
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
					for (int k=0; k<nn; ++k)
					{
					FENode& node = m.Node(fn[k]);
						
						if (fabs(node.m_r0.x) == m_bb_x) m_BN[fn[k]] = 1;
						if (fabs(node.m_r0.y) == m_bb_y) m_BN[fn[k]] = 1;
						if (fabs(node.m_r0.z) == m_bb_z) m_BN[fn[k]] = 1;
					}
				}
			}
		}
	}

}

//-----------------------------------------------------------------------------
bool FEMicroMaterial2O::PrepDisplacementBC()
{
	FEMesh& m = m_rve.GetMesh();
	int N = m.Nodes();

	// count the nr of exterior nodes
	int NN = 0, i;
	for (i=0; i<N; ++i) if (m_BN[i] == 1) ++NN;
	m_num_ext_node = NN;

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
bool FEMicroMaterial2O::PrepPeriodicBC()
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
void FEMicroMaterial2O::UpdateBC(FEModel& rve, mat3d& F, tens3drs& G)
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
		
		// LTE - Apply the second order boundary conditions to the RVE problem
		vec3d r1 = F*r0 + G.contractdyad1(r0)*0.5;
		//vec3d r1 = F*r0;

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
			vec3d u1 = U*u0 + G.contractdyad1(r1)/2. - G.contractdyad1(r0)/2.;

			// set this as the scale parameter for the offset
			FEParam* pp = pc->GetParameterList().Find("offset");
			assert(pp);
			pp->m_vscl = u1;
			pp->m_nlc = 0;
		}
	}
}

//-----------------------------------------------------------------------------
// The stiffness is evaluated at the same time the stress is evaluated so we 
// can just return it here. Note that this assumes that the stress function 
// is always called prior to the tangent function.
tens4ds FEMicroMaterial2O::Tangent(FEMaterialPoint &mp)
{
	FEMicroMaterialPoint2O& mmpt = *mp.ExtractData<FEMicroMaterialPoint2O>();
	return mmpt.m_Ca;
}

//-----------------------------------------------------------------------------
mat3ds FEMicroMaterial2O::Stress(FEMaterialPoint &mp)
{
	//// get the deformation gradient
	//FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	//FEMicroMaterialPoint2O& pt2O = *mp.ExtractData<FEMicroMaterialPoint2O>();
	//mat3d F = pt.m_F;
	//tens3drs G = pt2O.m_G;

	//// Create a local copy of the rve
	//FEModel rve;
	//rve.CopyFrom(m_rve);
	//rve.GetStep(0)->SetPrintLevel(FE_PRINT_NEVER);

	//// initialize
	//if (rve.Init() == false) throw FEMultiScaleException();

	//// apply the BC's
	//UpdateBC(rve, F, G);

	//// solve the RVE
	//bool bret = rve.Solve();

	//// set the plot file
	//FEBioPlotFile* pplt = new FEBioPlotFile(rve);
	//vector<int> item;
	//pplt->AddVariable("displacement", item);
	//pplt->AddVariable("stress", item);
	//pplt->Open(rve, "rve.xplt");
	//pplt->Write(rve);
	//pplt->Close();

	//// make sure it converged
	//if (bret == false) throw FEMultiScaleException();

	//// calculate the averaged stress
	//mat3ds sa = AveragedStress(rve, mp);

	//// calculate the averaged stiffness
	//AveragedStiffness(rve, mp, pt2O.m_Ca, pt2O.m_Da, pt2O.m_Ea);
	
	mat3ds sa; sa.zero();
	return sa;
}

//-----------------------------------------------------------------------------
//! Calculate the average stress from the RVE solution.
mat3ds FEMicroMaterial2O::AveragedStress(FEModel& rve, FEMaterialPoint &mp)
{
	//FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	//mat3d F = pt.m_F;
	//double J = pt.m_J;

	//// get the RVE mesh
	//FEMesh& m = rve.GetMesh();

	//mat3d T; T.zero();

	//// for periodic BC's we take the reaction forces directly from the periodic constraints
	//if (m_bperiodic)
	//{
	//	// get the reaction for from the periodic constraints
	//	for (int i=0; i<3; ++i)
	//	{
	//		FEPeriodicBoundary* pbc = dynamic_cast<FEPeriodicBoundary*>(rve.SurfacePairInteraction(i));
	//		assert(pbc);
	//		FEPeriodicSurface& ss = pbc->m_ss;
	//		int N = ss.Nodes();
	//		for (int i=0; i<N; ++i)
	//		{
	//			FENode& node = ss.Node(i);
	//			vec3d f = ss.m_Fr[i];

	//			// We multiply by two since the reaction forces are only stored at the slave surface 
	//			// and we also need to sum over the master nodes (NOTE: should I figure out a way to 
	//			// store the reaction forces on the master nodes as well?)
	//			T += (f & node.m_rt)*2.0;
	//		}
	//	}
	//}

	//// get the reaction force vector from the solid solver
	//// (We also need to do this for the periodic BC, since at the prescribed nodes,
	//// the contact forces will be zero). 
	//FEAnalysis* pstep = rve.GetCurrentStep();
	//FESolidSolver* ps = dynamic_cast<FESolidSolver*>(pstep->m_psolver);
	//assert(ps);
	//vector<double>& R = ps->m_Fr;
	//int nbc = rve.PrescribedBCs();
	//for (int i=0; i<nbc/3; ++i)
	//{
	//	FEPrescribedBC& dc = *rve.PrescribedBC(3*i);
	//	FENode& n = m.Node(dc.node);
	//	vec3d f;
	//	f.x = R[-n.m_ID[DOF_X]-2];
	//	f.y = R[-n.m_ID[DOF_Y]-2];
	//	f.z = R[-n.m_ID[DOF_Z]-2];
	//	T += f & n.m_rt;
	//}

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

	// Create a local copy of the rve
	FEModel rve;
	rve.CopyFrom(m_rve);
	rve.GetStep(0)->SetPrintLevel(FE_PRINT_NEVER);

	// initialize
	if (rve.Init() == false) throw FEMultiScaleException();

	// apply the BC's
	UpdateBC(rve, F, G);

	// solve the RVE
	bool bret = rve.Solve();

	// set the plot file
	FEBioPlotFile* pplt = new FEBioPlotFile(rve);
	vector<int> item;
	pplt->AddVariable("displacement", item);
	pplt->AddVariable("stress", item);
	pplt->Open(rve, "rve.xplt");
	pplt->Write(rve);
	pplt->Close();

	// make sure it converged
	if (bret == false) throw FEMultiScaleException();

	// calculate the averaged stress
	mat3ds sa; sa.zero();
	tens3ds taua; taua.zero();

	AveragedStress2O(rve, mp, pt.m_s, mmpt2O.m_tau);
	AveragedStress2OPK1(rve, mp, mmpt2O.m_PK1, mmpt2O.m_QK1);
	AveragedStress2OPK2(rve, mp, mmpt2O.m_S, mmpt2O.m_T);

	// calculate the averaged stiffness
	AveragedStiffness(rve, mp, mmpt2O.m_Ca, mmpt2O.m_Da, mmpt2O.m_Ea);

	calc_energy_diff(rve, mp, sa, taua);	
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

	mat3d T; T.zero();
	tens3ds tau; tau.zero();

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

				vec3d x; x = node.m_rt;
		
				tau.d[0] += (x.x*f.x*x.x)*2.0; 
				tau.d[1] += ((x.x*f.x*x.y + x.x*f.y*x.x + x.y*f.x*x.x)/3.)*2.0; 
				tau.d[2] += ((x.x*f.x*x.z + x.x*f.z*x.x + x.z*f.x*x.x)/3.)*2.0;
				tau.d[3] += ((x.x*f.y*x.y + x.y*f.x*x.y + x.y*f.y*x.x)/3.)*2.0; 
				tau.d[4] += ((x.x*f.y*x.z + x.y*f.x*x.z + x.z*f.y*x.x + x.x*f.z*x.y + x.z*f.x*x.y + x.y*f.z*x.x)/6.)*2.0; 
				tau.d[5] += ((x.x*f.z*x.z + x.z*f.x*x.z + x.z*f.z*x.x)/3.)*2.0;
				tau.d[6] += (x.y*f.y*x.y)*2.0; 
				tau.d[7] += ((x.y*f.y*x.z + x.y*f.z*x.y + x.z*f.y*x.y)/3.)*2.0;
				tau.d[8] += ((x.y*f.z*x.z + x.z*f.y*x.z + x.z*f.z*x.y)/3.)*2.0;
				tau.d[9] += (x.z*f.z*x.z)*2.0;
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

		vec3d x; x = n.m_rt;
		
		tau.d[0] += x.x*f.x*x.x; 
		tau.d[1] += (x.x*f.x*x.y + x.x*f.y*x.x + x.y*f.x*x.x)/3.; 
		tau.d[2] += (x.x*f.x*x.z + x.x*f.z*x.x + x.z*f.x*x.x)/3.;
		tau.d[3] += (x.x*f.y*x.y + x.y*f.x*x.y + x.y*f.y*x.x)/3.; 
		tau.d[4] += (x.x*f.y*x.z + x.y*f.x*x.z + x.z*f.y*x.x + x.x*f.z*x.y + x.z*f.x*x.y + x.y*f.z*x.x)/6.; 
		tau.d[5] += (x.x*f.z*x.z + x.z*f.x*x.z + x.z*f.z*x.x)/3.;
		tau.d[6] += x.y*f.y*x.y; 
		tau.d[7] += (x.y*f.y*x.z + x.y*f.z*x.y + x.z*f.y*x.y)/3.;
		tau.d[8] += (x.y*f.z*x.z + x.z*f.y*x.z + x.z*f.z*x.y)/3.;
		tau.d[9] += x.z*f.z*x.z; 
	}

	sa = T.sym() / (J*m_V0);
	taua = tau / 2*(J*m_V0);
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
					if ((m_BN[el.m_node[i]] == 1) && (m_BN[el.m_node[j]] == 1))
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


//-----------------------------------------------------------------------------
void FEMicroMaterial2O::calculate_d2O(tens5ds& d, double K[3][3], double Ri[3], double Rj[3] )
{
	d.d[0] += 0.5*(Ri[0]*K[0][0]*Rj[0]*Rj[0] + Ri[0]*Ri[0]*K[0][0]*Rj[0]);
	d.d[1] += 0.5*(Ri[0]*K[0][0]*Rj[0]*Rj[1] + Ri[0]*Ri[0]*K[0][0]*Rj[1]);
	d.d[2] += 0.5*(Ri[0]*K[0][0]*Rj[0]*Rj[2] + Ri[0]*Ri[0]*K[0][0]*Rj[2]);
	d.d[3] += 0.5*(Ri[0]*K[0][0]*Rj[1]*Rj[1] + Ri[0]*Ri[0]*K[0][1]*Rj[1]);
	d.d[4] += 0.5*(Ri[0]*K[0][0]*Rj[1]*Rj[2] + Ri[0]*Ri[0]*K[0][1]*Rj[2]);
	d.d[5] += 0.5*(Ri[0]*K[0][0]*Rj[2]*Rj[2] + Ri[0]*Ri[0]*K[0][2]*Rj[2]);
	d.d[6] += 0.5*(Ri[0]*K[0][1]*Rj[1]*Rj[1] + Ri[0]*Ri[0]*K[1][1]*Rj[1]);
	d.d[7] += 0.5*(Ri[0]*K[0][1]*Rj[1]*Rj[2] + Ri[0]*Ri[0]*K[1][1]*Rj[2]);
	d.d[8] += 0.5*(Ri[0]*K[0][1]*Rj[2]*Rj[2] + Ri[0]*Ri[0]*K[1][2]*Rj[2]);
	d.d[9] += 0.5*(Ri[0]*K[0][2]*Rj[2]*Rj[2] + Ri[0]*Ri[0]*K[2][2]*Rj[2]);
	
	d.d[10] += 0.5*(Ri[0]*K[1][0]*Rj[1]*Rj[1] + Ri[0]*Ri[1]*K[0][1]*Rj[1]);
	d.d[11] += 0.5*(Ri[0]*K[1][0]*Rj[1]*Rj[2] + Ri[0]*Ri[1]*K[0][1]*Rj[2]);
	d.d[12] += 0.5*(Ri[0]*K[1][0]*Rj[2]*Rj[2] + Ri[0]*Ri[1]*K[0][2]*Rj[2]);
	d.d[13] += 0.5*(Ri[0]*K[1][1]*Rj[1]*Rj[1] + Ri[0]*Ri[1]*K[1][1]*Rj[1]);
	d.d[14] += 0.5*(Ri[0]*K[1][1]*Rj[1]*Rj[2] + Ri[0]*Ri[1]*K[1][1]*Rj[2]);
	d.d[15] += 0.5*(Ri[0]*K[1][1]*Rj[2]*Rj[2] + Ri[0]*Ri[1]*K[1][2]*Rj[2]);
	d.d[16] += 0.5*(Ri[0]*K[1][2]*Rj[2]*Rj[2] + Ri[0]*Ri[1]*K[2][2]*Rj[2]);

	d.d[17] += 0.5*(Ri[0]*K[2][0]*Rj[2]*Rj[2] + Ri[0]*Ri[2]*K[0][2]*Rj[2]);
	d.d[18] += 0.5*(Ri[0]*K[2][1]*Rj[1]*Rj[1] + Ri[0]*Ri[2]*K[1][1]*Rj[1]);
	d.d[19] += 0.5*(Ri[0]*K[2][1]*Rj[1]*Rj[2] + Ri[0]*Ri[2]*K[1][1]*Rj[2]);
	d.d[20] += 0.5*(Ri[0]*K[2][1]*Rj[2]*Rj[2] + Ri[0]*Ri[2]*K[1][2]*Rj[2]);
	d.d[21] += 0.5*(Ri[0]*K[2][2]*Rj[2]*Rj[2] + Ri[0]*Ri[2]*K[2][2]*Rj[2]);

	d.d[22] += 0.5*(Ri[1]*K[1][1]*Rj[1]*Rj[1] + Ri[1]*Ri[1]*K[1][1]*Rj[1]);
	d.d[23] += 0.5*(Ri[1]*K[1][1]*Rj[1]*Rj[2] + Ri[1]*Ri[1]*K[1][1]*Rj[2]);
	d.d[24] += 0.5*(Ri[1]*K[1][1]*Rj[2]*Rj[2] + Ri[1]*Ri[1]*K[1][2]*Rj[2]);
	d.d[25] += 0.5*(Ri[1]*K[1][2]*Rj[2]*Rj[2] + Ri[1]*Ri[1]*K[2][2]*Rj[2]);

	d.d[26] += 0.5*(Ri[1]*K[2][1]*Rj[2]*Rj[2] + Ri[1]*Ri[2]*K[1][2]*Rj[2]);

	d.d[27] += 0.5*(Ri[1]*K[2][2]*Rj[2]*Rj[2] + Ri[1]*Ri[2]*K[2][2]*Rj[2]);

	d.d[28] += 0.5*(Ri[2]*K[2][2]*Rj[2]*Rj[2] + Ri[2]*Ri[2]*K[2][2]*Rj[2]);
}

//-----------------------------------------------------------------------------
void FEMicroMaterial2O::calculate_e2O(tens6ds& e, double K[3][3], double Ri[3], double Rj[3] )
{
	e.d[0] += Ri[0]*Ri[0]*K[0][0]*Rj[0]*Rj[0];
	e.d[1] += Ri[0]*Ri[0]*K[0][0]*Rj[0]*Rj[1];
	e.d[2] += Ri[0]*Ri[0]*K[0][0]*Rj[0]*Rj[2];
	e.d[3] += Ri[0]*Ri[0]*K[0][0]*Rj[1]*Rj[1];
	e.d[4] += Ri[0]*Ri[0]*K[0][0]*Rj[1]*Rj[2];
	e.d[5] += Ri[0]*Ri[0]*K[0][0]*Rj[2]*Rj[2];
	e.d[6] += Ri[0]*Ri[0]*K[0][1]*Rj[1]*Rj[1];
	e.d[7] += Ri[0]*Ri[0]*K[0][1]*Rj[1]*Rj[2];
	e.d[8] += Ri[0]*Ri[0]*K[0][1]*Rj[2]*Rj[2];
	e.d[9] += Ri[0]*Ri[0]*K[0][2]*Rj[2]*Rj[2];
	
	e.d[10] += Ri[0]*Ri[0]*K[1][0]*Rj[1]*Rj[1];
	e.d[11] += Ri[0]*Ri[0]*K[1][0]*Rj[1]*Rj[2];
	e.d[12] += Ri[0]*Ri[0]*K[1][0]*Rj[2]*Rj[2];
	e.d[13] += Ri[0]*Ri[0]*K[1][1]*Rj[1]*Rj[1];
	e.d[14] += Ri[0]*Ri[0]*K[1][1]*Rj[1]*Rj[2];
	e.d[15] += Ri[0]*Ri[0]*K[1][1]*Rj[2]*Rj[2];
	e.d[16] += Ri[0]*Ri[0]*K[1][2]*Rj[2]*Rj[2];

	e.d[17] += Ri[0]*Ri[0]*K[2][0]*Rj[1]*Rj[1];
	e.d[18] += Ri[0]*Ri[0]*K[2][0]*Rj[1]*Rj[2];
	e.d[19] += Ri[0]*Ri[0]*K[2][0]*Rj[2]*Rj[2];
	e.d[20] += Ri[0]*Ri[0]*K[2][1]*Rj[1]*Rj[1];
	e.d[21] += Ri[0]*Ri[0]*K[2][1]*Rj[1]*Rj[2];
	e.d[22] += Ri[0]*Ri[0]*K[2][1]*Rj[2]*Rj[2];
	e.d[23] += Ri[0]*Ri[0]*K[2][2]*Rj[2]*Rj[2];

	e.d[24] += Ri[0]*Ri[1]*K[1][0]*Rj[2]*Rj[2];
	e.d[25] += Ri[0]*Ri[1]*K[1][1]*Rj[1]*Rj[1];
	e.d[26] += Ri[0]*Ri[1]*K[1][1]*Rj[1]*Rj[2];
	e.d[27] += Ri[0]*Ri[1]*K[1][1]*Rj[2]*Rj[2];
	e.d[28] += Ri[0]*Ri[1]*K[1][2]*Rj[2]*Rj[2];

	e.d[29] += Ri[0]*Ri[1]*K[2][0]*Rj[2]*Rj[2];
	e.d[30] += Ri[0]*Ri[1]*K[2][1]*Rj[1]*Rj[1];
	e.d[31] += Ri[0]*Ri[1]*K[2][1]*Rj[1]*Rj[2];
	e.d[32] += Ri[0]*Ri[1]*K[2][1]*Rj[2]*Rj[2];
	e.d[33] += Ri[0]*Ri[1]*K[2][2]*Rj[2]*Rj[2];

	e.d[34] += Ri[0]*Ri[2]*K[2][1]*Rj[1]*Rj[1];
	e.d[35] += Ri[0]*Ri[2]*K[2][1]*Rj[1]*Rj[2];
	e.d[36] += Ri[0]*Ri[2]*K[2][1]*Rj[2]*Rj[2];
	e.d[37] += Ri[0]*Ri[2]*K[2][2]*Rj[2]*Rj[2];

	e.d[38] += Ri[1]*Ri[1]*K[1][1]*Rj[1]*Rj[1];
	e.d[39] += Ri[1]*Ri[1]*K[1][1]*Rj[1]*Rj[2];
	e.d[40] += Ri[1]*Ri[1]*K[1][1]*Rj[2]*Rj[2];
	e.d[41] += Ri[1]*Ri[1]*K[1][2]*Rj[2]*Rj[2];

	e.d[42] += Ri[1]*Ri[1]*K[2][1]*Rj[2]*Rj[2];
	e.d[43] += Ri[1]*Ri[1]*K[2][2]*Rj[2]*Rj[2];

	e.d[44] += Ri[1]*Ri[2]*K[2][2]*Rj[2]*Rj[2];

	e.d[45] += Ri[2]*Ri[2]*K[2][2]*Rj[2]*Rj[2];
}


//-----------------------------------------------------------------------------
//! Calculate the energy difference between the RVE problem and the macro material point
void FEMicroMaterial2O::calc_energy_diff(FEModel& rve, FEMaterialPoint& mp, mat3ds& sa, tens3ds& taua)
{
	// get the deformation gradient and deformation hessian
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	FEMicroMaterialPoint2O& mmpt2O = *mp.ExtractData<FEMicroMaterialPoint2O>();
	mat3d F = pt.m_F;
	tens3drs G = mmpt2O.m_G;

	mat3d Ftrans = F.transpose();
	tens3dls Gtrans = G.transpose();

	mat3d Finv = F.inverse();
	mat3d Finvtrans = Finv.transpose();
	tens3drs Ginv; Ginv = G; Ginv.contractleg2(Finv,1); Ginv.contractleg2(Finv,2); Ginv.contractleg2(Finv,3);
	tens3dls Ginvtrans = Ginv.transpose();
	
	// calculate infinitesimal strain
	mmpt2O.m_inf_str = ((F.transpose() + F)*0.5 - mat3dd(1)).sym();
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
	mmpt2O.m_H = ((Gtrans.multiply2right(F).LStoUnsym() + G.multiply2left(Ftrans).RStoUnsym())*0.5).symm();

	// calculate Euler-Almansi strain
	mmpt2O.m_e = ((mat3dd(1) - Finvtrans*Finv)*0.5).sym();
	mmpt2O.m_h = ((Ginvtrans.multiply2right(Finv).LStoUnsym() + Ginv.multiply2left(Finvtrans).RStoUnsym())*-0.5).symm();
	
	// calculate the energy difference between macro point and RVE
	// to verify that we have satisfied the Hill-Mandel condition
	double macro_energy = sa.dotdot(mmpt2O.m_e) + taua.tripledot3s(mmpt2O.m_h);
	
	double rve_energy_avg = 0.;
	int nint; 
	double* w;
	double v = 0.;

	FEMesh& m = rve.GetMesh();
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
				FEElasticMaterialPoint& rve_pt = *el.GetMaterialPoint(n)->ExtractData<FEElasticMaterialPoint>();
				mat3d rve_F = rve_pt.m_F;
				mat3ds rve_e = ((mat3dd(1) - rve_F.transinv()*rve_F.inverse())*0.5).sym();
				mat3ds rve_s = rve_pt.m_s;
				rve_energy_avg += rve_s.dotdot(rve_e)*rve_pt.m_J*w[n];
				v += rve_pt.m_J*w[n];
			}
		}
	}

	rve_energy_avg /= v;
	mmpt2O.m_energy_diff = fabs(macro_energy - rve_energy_avg);
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
				PK1 += (f & node.m_r0)*2.0;

				vec3d X; X = node.m_r0;
		
				QK1.d[0] += 2.0*f.x*X.x*X.x; 
				QK1.d[1] += 2.0*f.x*X.x*X.y;
				QK1.d[2] += 2.0*f.x*X.x*X.z;
				QK1.d[3] += 2.0*f.x*X.y*X.y;
				QK1.d[4] += 2.0*f.x*X.y*X.z;
				QK1.d[5] += 2.0*f.x*X.z*X.z;
				QK1.d[6] += 2.0*f.y*X.x*X.x;
				QK1.d[7] += 2.0*f.y*X.x*X.y;
				QK1.d[8] += 2.0*f.y*X.x*X.z;
				QK1.d[9] += 2.0*f.y*X.y*X.y;
				QK1.d[10] += 2.0*f.y*X.y*X.z;
				QK1.d[11] += 2.0*f.y*X.z*X.z;
				QK1.d[12] += 2.0*f.z*X.x*X.x;
				QK1.d[13] += 2.0*f.z*X.x*X.y;
				QK1.d[14] += 2.0*f.z*X.x*X.z;
				QK1.d[15] += 2.0*f.z*X.y*X.y;
				QK1.d[16] += 2.0*f.z*X.y*X.z;
				QK1.d[17] += 2.0*f.z*X.z*X.z;
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
		
		PK1 += f & n.m_r0;
		vec3d X; X = n.m_r0;
		
		QK1.d[0] += f.x*X.x*X.x; 
		QK1.d[1] += f.x*X.x*X.y;
		QK1.d[2] += f.x*X.x*X.z;
		QK1.d[3] += f.x*X.y*X.y;
		QK1.d[4] += f.x*X.y*X.z;
		QK1.d[5] += f.x*X.z*X.z;
		QK1.d[6] += f.y*X.x*X.x;
		QK1.d[7] += f.y*X.x*X.y;
		QK1.d[8] += f.y*X.x*X.z;
		QK1.d[9] += f.y*X.y*X.y;
		QK1.d[10] += f.y*X.y*X.z;
		QK1.d[11] += f.y*X.z*X.z;
		QK1.d[12] += f.z*X.x*X.x;
		QK1.d[13] += f.z*X.x*X.y;
		QK1.d[14] += f.z*X.x*X.z;
		QK1.d[15] += f.z*X.y*X.y;
		QK1.d[16] += f.z*X.y*X.z;
		QK1.d[17] += f.z*X.z*X.z;
	}

	PK1a = PK1 / m_V0;
	QK1a = QK1 / (2*m_V0);
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
			FEPeriodicBoundary* pbc = dynamic_cast<FEPeriodicBoundary*>(rve.SurfacePairInteraction(i));
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

				vec3d X; X = node.m_r0;
		
				T.d[0] += (X.x*f0.x*X.x)*2.0; 
				T.d[1] += ((X.x*f0.x*X.y + X.x*f0.y*X.x + X.y*f0.x*X.x)/3.)*2.0; 
				T.d[2] += ((X.x*f0.x*X.z + X.x*f0.z*X.x + X.z*f0.x*X.x)/3.)*2.0;
				T.d[3] += ((X.x*f0.y*X.y + X.y*f0.x*X.y + X.y*f0.y*X.x)/3.)*2.0; 
				T.d[4] += ((X.x*f0.y*X.z + X.y*f0.x*X.z + X.z*f0.y*X.x + X.x*f0.z*X.y + X.z*f0.x*X.y + X.y*f0.z*X.x)/6.)*2.0; 
				T.d[5] += ((X.x*f0.z*X.z + X.z*f0.x*X.z + X.z*f0.z*X.x)/3.)*2.0;
				T.d[6] += (X.y*f0.y*X.y)*2.0; 
				T.d[7] += ((X.y*f0.y*X.z + X.y*f0.z*X.y + X.z*f0.y*X.y)/3.)*2.0;
				T.d[8] += ((X.y*f0.z*X.z + X.z*f0.y*X.z + X.z*f0.z*X.y)/3.)*2.0;
				T.d[9] += (X.z*f0.z*X.z)*2.0;
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
		vec3d f0 = Finv*f;
		
		S += f0 & n.m_r0;

		vec3d X; X = n.m_r0;
		
		T.d[0] += X.x*f0.x*X.x; 
		T.d[1] += (X.x*f0.x*X.y + X.x*f0.y*X.x + X.y*f0.x*X.x)/3.; 
		T.d[2] += (X.x*f0.x*X.z + X.x*f0.z*X.x + X.z*f0.x*X.x)/3.;
		T.d[3] += (X.x*f0.y*X.y + X.y*f0.x*X.y + X.y*f0.y*X.x)/3.; 
		T.d[4] += (X.x*f0.y*X.z + X.y*f0.x*X.z + X.z*f0.y*X.x + X.x*f0.z*X.y + X.z*f0.x*X.y + X.y*f0.z*X.x)/6.; 
		T.d[5] += (X.x*f0.z*X.z + X.z*f0.x*X.z + X.z*f0.z*X.x)/3.;
		T.d[6] += X.y*f0.y*X.y; 
		T.d[7] += (X.y*f0.y*X.z + X.y*f0.z*X.y + X.z*f0.y*X.y)/3.;
		T.d[8] += (X.y*f0.z*X.z + X.z*f0.y*X.z + X.z*f0.z*X.y)/3.;
		T.d[9] += X.z*f0.z*X.z; 
	}

	Sa = S.sym() / m_V0;
	Ta = T / (2*m_V0);
}