/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in 
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
#include "FEMicroMaterial.h"
#include "FECore/FEElemElemList.h"
#include "FECore/log.h"
#include "FESolidSolver2.h"
#include "FEElasticSolidDomain.h"
#include "FECore/FEAnalysis.h"
//#include "FEBioXML/FEBioImport.h"
//#include "FEBioPlot/FEBioPlotFile.h"
#include <FECore/mat6d.h>
#include "FEBCPrescribedDeformation.h"
#include <sstream>

//=============================================================================
FERVEProbe::FERVEProbe(FEModel& fem, FEModel& rve, const char* szfile) : FECallBack(&fem, CB_ALWAYS), m_rve(rve), m_file(szfile) 
{
	m_bdebug = false;
}

bool FERVEProbe::Execute(FEModel& fem, int nwhen)
{
	if (nwhen == CB_INIT)	// initialize the plot file
	{
		// create a plot file
/*		m_xplt = new FEBioPlotFile(m_rve);
		if (m_xplt->Open(m_rve, m_file.c_str()) == false)
		{
			feLog("Failed creating probe.\n\n");
			delete m_xplt; m_xplt = 0;
		}
*/
		// write the initial state
		Save();
	}
	else if (nwhen == CB_MINOR_ITERS)
	{
		if (m_bdebug) Save();
	}
	else if (nwhen == CB_MAJOR_ITERS)	// store the current state
	{
		Save();
	}
	else if (nwhen == CB_SOLVED)	// clean up
	{
		if (m_xplt) delete m_xplt;
		m_xplt = 0;
	}

	return true;
}

void FERVEProbe::Save()
{
//	if (m_xplt) m_xplt->Write(m_rve, (float) m_rve.GetCurrentTime());
}

//=============================================================================
FEMicroMaterialPoint::FEMicroMaterialPoint(FEMaterialPoint* mp) : FEMaterialPoint(mp)
{
	m_macro_energy = 0.;
	m_micro_energy = 0.;
	m_energy_diff = 0.;
	
	m_macro_energy_inc = 0.;
	m_micro_energy_inc = 0.;
}

//-----------------------------------------------------------------------------
//! Initialize material point data
void FEMicroMaterialPoint::Init()
{
	FEMaterialPoint::Init();
	m_F_prev.unit();
}

//-----------------------------------------------------------------------------
//! Initialize material point data
void FEMicroMaterialPoint::Update(const FETimeInfo& timeInfo)
{
	FEMaterialPoint::Update(timeInfo);
	FEElasticMaterialPoint& pt = *ExtractData<FEElasticMaterialPoint>();
	m_F_prev = pt.m_F;
}

//-----------------------------------------------------------------------------
//! create a shallow copy
FEMaterialPoint* FEMicroMaterialPoint::Copy()
{
	FEMicroMaterialPoint* pt = new FEMicroMaterialPoint(m_pNext?m_pNext->Copy():0);
	return pt;
}

//-----------------------------------------------------------------------------
//! serialize material point data
void FEMicroMaterialPoint::Serialize(DumpStream& ar)
{
	FEMaterialPoint::Serialize(ar);
	if (ar.IsSaving())
	{
	}
	else
	{
	}
}

//=============================================================================
BEGIN_FECORE_CLASS(FEMicroProbe, FEMaterial)
	ADD_PARAMETER(m_neid  , "element_id");
	ADD_PARAMETER(m_ngp   , "gausspt"   );
	ADD_PARAMETER(m_szfile, "file"      );
	ADD_PARAMETER(m_bdebug, "debug"     );
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEMicroProbe::FEMicroProbe(FEModel* pfem) : FEMaterial(pfem)
{
	m_neid = -1;	// invalid element - this must be defined by user
	m_ngp = 1;		// by default, first gauss point (note is one-based!)
	m_szfile = "rve.xplt";
	m_probe = 0;
	m_bdebug = false;
}

FEMicroProbe::~FEMicroProbe()
{
	if (m_probe) delete m_probe;
}

//=============================================================================

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEMicroMaterial, FEElasticMaterial)
	ADD_PARAMETER(m_szrve    , "RVE"     );
	ADD_PARAMETER(m_szbc     , "bc_set"  );
	ADD_PARAMETER(m_bctype   , "rve_type" );
	ADD_PARAMETER(m_scale	 , "scale"   ); 

	ADD_PROPERTY(m_probe, "probe", false);

END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEMicroMaterial::FEMicroMaterial(FEModel* pfem) : FEElasticMaterial(pfem)
{
	// initialize parameters
	m_szrve[0] = 0;
	m_szbc[0] = 0;
	m_bctype = FERVEModel::DISPLACEMENT;	// use displacement BCs by default
	m_scale = 1.0;
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
bool FEMicroMaterial::Init()
{
	if (FEElasticMaterial::Init() == false) return false;

	// load the RVE model
/*	FEBioImport fim;
	if (fim.Load(m_mrve, m_szrve.c_str()) == false)
	{
		return false;
	}
*/
	// We don't want to output anything from the RVE
	m_mrve.BlockLog();

	// scale the RVE
	if (m_scale != 1.0) m_mrve.ScaleGeometry(m_scale);

	// initialize the RVE model
	// This also creates the necessary boundary conditions
	bool bret = m_mrve.InitRVE(m_bctype, m_szbc.c_str()); 

	if (bret == false) {
		feLogError("An error occurred preparing RVE model"); return false;
	}

	return true;
}

//-----------------------------------------------------------------------------
// Note that this function is not used in the first-order implemenetation
mat3ds FEMicroMaterial::Stress(FEMaterialPoint &mp)
{
	// get the deformation gradient
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	FEMicroMaterialPoint& mmpt = *mp.ExtractData<FEMicroMaterialPoint>();
	mat3d F = pt.m_F;
	
	// update the BC's
	mmpt.m_rve.Update(F);

	// solve the RVE
	bool bret = mmpt.m_rve.Solve();

	// make sure it converged
	if (bret == false) throw FEMultiScaleException(-1, -1);

	// calculate the averaged Cauchy stress
	mat3ds sa = mmpt.m_rve.StressAverage(mp);
	
	// calculate the difference between the macro and micro energy for Hill-Mandel condition
	mmpt.m_micro_energy = micro_energy(mmpt.m_rve);	
	
	return sa;
}

//-----------------------------------------------------------------------------
// The stiffness is evaluated at the same time the stress is evaluated so we 
// can just return it here. Note that this assumes that the stress function 
// is always called prior to the tangent function.
tens4ds FEMicroMaterial::Tangent(FEMaterialPoint &mp)
{
	FEMicroMaterialPoint& mmpt = *mp.ExtractData<FEMicroMaterialPoint>();
	return mmpt.m_rve.StiffnessAverage(mp);
}

//-----------------------------------------------------------------------------
//! Calculate the "energy" of the RVE model, i.e. the volume averaged of PK1:F
double FEMicroMaterial::micro_energy(FEModel& rve)
{
	double E_avg = 0.0;
	double V0 = 0.0;
	FEMesh& m = rve.GetMesh();
	for (int k=0; k<m.Domains(); ++k)
	{
		FESolidDomain& dom = static_cast<FESolidDomain&>(m.Domain(k));
		for (int i=0; i<dom.Elements(); ++i)
		{
			FESolidElement& el = dom.Element(i);
			int nint = el.GaussPoints();
			double* w = el.GaussWeights();
			
			for (int n=0; n<nint; ++n)
			{
				FEMaterialPoint& mp = *el.GetMaterialPoint(n);
				FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

				mat3d& F = pt.m_F;
				double J = F.det();
				mat3ds& s = pt.m_s;	// Cauchy stress

				// PK1 stress
				mat3d Pk1 = J*s*F.transinv();

				double energy = Pk1.dotdot(F);

				double J0 = dom.detJ0(el, n);		
				E_avg += energy*w[n]*J0;

				V0 += w[n]*J0;
			}
		}
	}

	return E_avg/V0;
}

//-----------------------------------------------------------------------------
//! Calculate the average stress from the RVE solution.
mat3d FEMicroMaterial::AveragedStressPK1(FEModel& rve, FEMaterialPoint &mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	mat3d F = pt.m_F;
	double J = pt.m_J;
	
	// get the RVE mesh
	FEMesh& m = rve.GetMesh();

	mat3d PK1; PK1.zero();

	// for periodic BC's we take the reaction forces directly from the periodic constraints
	if (m_bctype == FERVEModel::PERIODIC_AL)
	{
		// get the reaction for from the periodic constraints
		for (int i=0; i<3; ++i)
		{
			FEPeriodicBoundary1O* pbc = dynamic_cast<FEPeriodicBoundary1O*>(rve.SurfacePairConstraint(i));
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
	FEBCPrescribedDeformation& dc = dynamic_cast<FEBCPrescribedDeformation&>(*rve.BoundaryCondition(0));

	const FENodeSet& nset = *dc.GetNodeSet();
	int nitems = nset.Size();
	for (int i=0; i<nitems; ++i)
	{
		const FENode& n = *nset.Node(i);
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
mat3ds FEMicroMaterial::AveragedStressPK2(FEModel& rve, FEMaterialPoint &mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	mat3d F = pt.m_F;
	double J = pt.m_J;
	mat3d Finv = F.inverse();

	// get the RVE mesh
	FEMesh& m = rve.GetMesh();

	mat3d S; S.zero();

	// for periodic BC's we take the reaction forces directly from the periodic constraints
	if (m_bctype == FERVEModel::PERIODIC_AL)
	{
		// get the reaction for from the periodic constraints
		for (int i=0; i<3; ++i)
		{
			FEPeriodicBoundary1O* pbc = dynamic_cast<FEPeriodicBoundary1O*>(rve.SurfacePairConstraint(i));
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
	FEBCPrescribedDeformation& dc = dynamic_cast<FEBCPrescribedDeformation&>(*rve.BoundaryCondition(0));
	const FENodeSet& nset = *dc.GetNodeSet();
	int nitems = nset.Size();
	for (int i=0; i<nitems; ++i)
	{
		const FENode& n = *nset.Node(i);
		vec3d f;
		f.x = R[-n.m_ID[dof_X]-2];
		f.y = R[-n.m_ID[dof_Y]-2];
		f.z = R[-n.m_ID[dof_Z]-2];
		vec3d f0 = Finv*f;
		S += f0 & n.m_r0;
	}

	double V0 = m_mrve.InitialVolume();
	return S.sym() / V0;
}
