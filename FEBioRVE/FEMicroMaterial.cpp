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
#include "FEMicroMaterial.h"
#include "FECore/FEElemElemList.h"
#include "FECore/log.h"
#include "FEBioMech/FESolidSolver2.h"
#include "FEBioMech/FEElasticSolidDomain.h"
#include "FECore/FEAnalysis.h"
#include "FEBioXML/FEBioImport.h"
#include <FECore/mat6d.h>
#include "FEBioMech/FEBCPrescribedDeformation.h"
#include "FERVEProbe.h"
#include <sstream>

//=============================================================================
FEMicroMaterialPoint::FEMicroMaterialPoint()
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
	FEElasticMaterialPoint::Init();
	m_F_prev.unit();
}

//-----------------------------------------------------------------------------
//! Initialize material point data
void FEMicroMaterialPoint::Update(const FETimeInfo& timeInfo)
{
	FEElasticMaterialPoint::Update(timeInfo);
	m_F_prev = m_F;

	// clear rewind stack so the next rewind won't overwrite current state
	m_rve.RCI_ClearRewindStack();
}

//-----------------------------------------------------------------------------
//! create a shallow copy
FEMaterialPointData* FEMicroMaterialPoint::Copy()
{
	FEMicroMaterialPoint* pt = new FEMicroMaterialPoint();
	if (m_pNext) pt->SetNext(m_pNext->Copy());
	return pt;
}

//-----------------------------------------------------------------------------
//! serialize material point data
void FEMicroMaterialPoint::Serialize(DumpStream& ar)
{
	FEMaterialPointData::Serialize(ar);
	ar & m_S & m_F_prev;
	ar & m_macro_energy;
	ar & m_micro_energy;
	ar & m_energy_diff;
	ar & m_macro_energy_inc;
	ar & m_micro_energy_inc;
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
FEMaterialPointData* FEMicroMaterial::CreateMaterialPointData()
{
	return new FEMicroMaterialPoint;
}

//-----------------------------------------------------------------------------
bool FEMicroMaterial::Init()
{
	if (FEElasticMaterial::Init() == false) return false;

	// load the RVE model
#ifndef FEBIOMECH_EXPORTS
	FEBioImport fim;
	if (fim.Load(m_mrve, m_szrve.c_str()) == false)
	{
		return false;
	}
#else
	return false;
#endif

	// set the parent FEM
	m_mrve.SetParentModel(GetFEModel());

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
	FEMicroMaterialPoint& pt = *mp.ExtractData<FEMicroMaterialPoint>();
	mat3d F = pt.m_F;

	// calculate the averaged Cauchy stress
	mat3ds sa = pt.m_rve.StressAverage(F, mp);
	
	// calculate the difference between the macro and micro energy for Hill-Mandel condition
	pt.m_micro_energy = micro_energy(pt.m_rve);	
	
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
				vec3d f = ss.m_data[i].m_Fr;

				// We multiply by two since the reaction forces are only stored at the primary surface 
				// and we also need to sum over the secondary nodes (NOTE: should I figure out a way to 
				// store the reaction forces on the secondary nodes as well?)
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
				vec3d f = ss.m_data[i].m_Fr;
				vec3d f0 = Finv*f;

				// We multiply by two since the reaction forces are only stored at the primary surface 
				// and we also need to sum over the secondary nodes (NOTE: should I figure out a way to 
				// store the reaction forces on the secondary nodes as well?)
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
