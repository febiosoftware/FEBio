#include "StdAfx.h"
#include "FEMicroMaterial.h"
#include "FEElemElemList.h"
#include "Log.h"
#include "Console.h"

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
//	m_rve.m_pStep->SetPlotLevel(FE_PLOT_NEVER);

	// make sure we are using the same linear solver as the parent FEM
	m_rve.m_nsolver = PARDISO_SOLVER;

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
	for (int i=0; i<m.SolidElements(); ++i, ++M)
	{
		FESolidElement& el = m.SolidElement(i);
		nf = m.Faces(el);
		for (int j=0; j<nf; ++j)
		{
			if (EEL.Neighbor(M, j) == 0)
			{
				int nn = m.GetFace(el, j, fn);

				// mark all nodes
				tag[ el.m_node[ fn[0] ] ] = 1;
				tag[ el.m_node[ fn[1] ] ] = 1;
				tag[ el.m_node[ fn[2] ] ] = 1;
				if (nn == 4) tag[ el.m_node[ fn[3] ] ] = 1;
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

	// reset the logfile mode
	log.SetMode(nmode);
}

//-----------------------------------------------------------------------------
mat3ds FEMicroMaterial::Stress(FEMaterialPoint &mp)
{
	// the logfile is a shared resource between the master FEM and the RVE
	// in order not to corrupt the logfile we don't print anything for 
	// the RVE problem.
	Logfile& log = GetLogfile();
	Logfile::MODE nmode = log.GetMode();
	GetLogfile().SetMode(Logfile::NEVER);

	// reset the RVE
	m_rve.Reset();

	// get the deformation gradient
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	mat3d F = pt.F;

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
	return AveragedStress();
}

//-----------------------------------------------------------------------------

mat3ds FEMicroMaterial::AveragedStress()
{
	mat3ds s(0);
	FEMesh& m = m_rve.m_mesh;
	double V = 0, ve;
	int nint, n, i;
	double* w, J;
	for (i=0; i<m.SolidElements(); ++i)
	{
		FESolidElement& el = m.SolidElement(i);
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

	return (s / V);
}

//-----------------------------------------------------------------------------
tens4ds FEMicroMaterial::Tangent(FEMaterialPoint &pt)
{
	return dynamic_cast<FESolidMaterial*>(m_rve.GetMaterial(0))->Tangent(pt);
}
