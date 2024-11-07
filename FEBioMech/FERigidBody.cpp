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
#include "FERigidBody.h"
#include <FECore/FEMaterial.h>
#include <FECore/FESolidDomain.h>
#include <FECore/FEModel.h>
#include "RigidBC.h"
#include "FERigidMaterial.h"
#include "FERigidSolidDomain.h"
#include "FERigidShellDomain.h"
#include <FECore/log.h>

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FERigidBody, FECoreBase)
	ADD_PARAMETER(m_Fr.x, "Fx");
	ADD_PARAMETER(m_Fr.y, "Fy");
	ADD_PARAMETER(m_Fr.z, "Fz");
	ADD_PARAMETER(m_Mr.x, "Mx");
	ADD_PARAMETER(m_Mr.y, "My");
	ADD_PARAMETER(m_Mr.z, "Mz");
	ADD_PARAMETER(m_euler, "euler");
	ADD_PARAMETER(m_r0, "initial_position");
	ADD_PARAMETER(m_rt, "position");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FERigidBody::FERigidBody(FEModel* pfem) : FECoreBase(pfem)
{
    m_bpofr = false;
	for (int i=0; i<6; ++i)
	{
		m_pDC[i] = 0;
		m_LM[i] = -1;
		m_BC[i] = DOF_OPEN;
	}
	m_prb = 0;

	// zero total displacements
	m_Ut[0] = m_Up[0] = 0;
	m_Ut[1] = m_Up[1] = 0;
	m_Ut[2] = m_Up[2] = 0;
	m_Ut[3] = m_Up[3] = 0;
	m_Ut[4] = m_Up[4] = 0;
	m_Ut[5] = m_Up[5] = 0;

    // initialize velocity and acceleration of center of mass
    m_vt = m_at = vec3d(0,0,0);
    
	// initialize orientation
	m_qt = quatd(0, vec3d(0,0,1));
	m_euler = vec3d(0,0,0);
    
    // initialize angular velocity and acceleration
    m_wt = m_alt = vec3d(0,0,0);

	// initialize reaction forces
	m_Fr = m_Fp = vec3d(0,0,0);
	m_Mr = m_Mp = vec3d(0,0,0);
}

//-----------------------------------------------------------------------------
FERigidBody::~FERigidBody()
{

}

//-----------------------------------------------------------------------------
//! Reset rigid body data (called from FEM::Reset)
void FERigidBody::Reset()
{
	// zero total displacements
	m_Ut[0] = m_Up[0] = 0;
	m_Ut[1] = m_Up[1] = 0;
	m_Ut[2] = m_Up[2] = 0;
	m_Ut[3] = m_Up[3] = 0;
	m_Ut[4] = m_Up[4] = 0;
	m_Ut[5] = m_Up[5] = 0;

    // initialize velocity and acceleration of center of mass
    m_vp = m_vt = vec3d(0,0,0);
    m_ap = m_at = vec3d(0,0,0);
    
	// initialize orientation
	m_qp = m_qt = quatd(0, vec3d(0,0,1));
	m_euler = vec3d(0,0,0);

    // initialize angular velocity and acceleration
    m_wp = m_wt = vec3d(0,0,0);
    m_alp = m_alt = vec3d(0,0,0);
    
    // initialize angular momentum and its time rate of change
    m_hp = m_ht = vec3d(0,0,0);
    m_dhp = m_dht = vec3d(0,0,0);
    
	// initialize center of mass
	m_rt = m_r0;

	// reset reaction force and torque
	m_Fr = vec3d(0,0,0);
	m_Mr = vec3d(0,0,0);
}

//-----------------------------------------------------------------------------
//! This function is called at the start of each time step and is used to update
//! some variables.
bool FERigidBody::Init()
{
	// clear reaction forces
	m_Fr = m_Mr = vec3d(0,0,0);

	// store previous state
	m_rp = m_rt;
    m_vp = m_vt;
    m_ap = m_at;
	m_qp = m_qt;
    m_wp = m_wt;
    m_alp = m_alt;
    m_hp = m_ht;
    m_dhp = m_dht;
	m_Up[0] = m_Ut[0];
	m_Up[1] = m_Ut[1];
	m_Up[2] = m_Ut[2];
	m_Up[3] = m_Ut[3];
	m_Up[4] = m_Ut[4];
	m_Up[5] = m_Ut[5];

	// zero incremental displacements
	m_du[0] = m_dul[0] = 0.0;
	m_du[1] = m_dul[1] = 0.0;
	m_du[2] = m_dul[2] = 0.0;
	m_du[3] = m_dul[3] = 0.0;
	m_du[4] = m_dul[4] = 0.0;
	m_du[5] = m_dul[5] = 0.0;

	return true;
}

//-----------------------------------------------------------------------------
//! Set the rigid body's center of mass directly
void FERigidBody::SetCOM(vec3d rc)
{
	m_r0 = m_rp = m_rt = rc;
}

//-----------------------------------------------------------------------------
//! Update the mass of the rigid body
void FERigidBody::UpdateMass()
{
	// total mass of rigid body
	m_mass = 0;

	// loop over all domains
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	for (int nd = 0; nd < mesh.Domains(); ++nd)
	{
		FEDomain& dom = mesh.Domain(nd);
		FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(dom.GetMaterial());
		if (pm && (pm->GetRigidBodyID() == m_nID))
		{
			FERigidSolidDomain* pbd = dynamic_cast<FERigidSolidDomain*>(&dom);
			if (pbd) m_mass += pbd->CalculateMass();

			FERigidShellDomain* psd = dynamic_cast<FERigidShellDomain*>(&dom);
			if (psd) m_mass += psd->CalculateMass();
		}
	}

	if (m_mass == 0)
	{
		feLogWarning("Zero mass for rigid body \"%s.\"", GetName().c_str());
	}
}

//-----------------------------------------------------------------------------
//! Calculates the rigid body's total mass, center of mass, and mass moment
//! of inertia about the center of mass
//!
void FERigidBody::UpdateCOM()
{
	// center of mass
	vec3d rc(0,0,0);

	// loop over all domains
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	for (int nd=0; nd < mesh.Domains(); ++nd)
	{
		FEDomain& dom = mesh.Domain(nd);
		FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(dom.GetMaterial());
		if (pm && (pm->GetRigidBodyID() == m_nID))
		{
			FERigidSolidDomain* pbd = dynamic_cast<FERigidSolidDomain*>(&dom);
			if (pbd) rc += pbd->CalculateCOM();

			FERigidShellDomain* psd = dynamic_cast<FERigidShellDomain*>(&dom);
			if (psd) rc += psd->CalculateCOM();
		}
	}

	// normalize com
	if (m_mass != 0) rc /= m_mass;

	// store com
	SetCOM(rc);
}

//-----------------------------------------------------------------------------
void FERigidBody::UpdateMOI()
{
	// initialize some data
	mat3d moi(0, 0, 0, 0, 0, 0, 0, 0, 0);    // mass moment of inertia about origin

	// loop over all domains
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	for (int nd = 0; nd < mesh.Domains(); ++nd)
	{
		FEDomain& dom = mesh.Domain(nd);
		FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(dom.GetMaterial());
		if (pm && (pm->GetRigidBodyID() == m_nID))
		{
			FERigidSolidDomain* pbd = dynamic_cast<FERigidSolidDomain*>(&dom);
			if (pbd) moi += pbd->CalculateMOI();

			FERigidShellDomain* psd = dynamic_cast<FERigidShellDomain*>(&dom);
			if (psd) moi += psd->CalculateMOI();
		}
	}

	// use parallel axis theorem to transfer moi to com
	// and store moi
	mat3dd I(1);        // identity tensor
	vec3d rc = m_r0;
	m_moi = moi.sym() - m_mass*((rc*rc)*I - dyad(rc));
}

//-----------------------------------------------------------------------------
vec3d FERigidBody::CayleyIncrementalCompoundRotation()
{
    // incremental rotation in spatial frame
    quatd q = m_qt*m_qp.Inverse();
    q.MakeUnit();                           // clean-up roundoff errors
    double theta = 2*tan(q.GetAngle()/2);   // get theta from Cayley transform
    vec3d e = q.GetVector();
    
    return e*theta;
}

//-----------------------------------------------------------------------------
void FERigidBody::Serialize(DumpStream& ar)
{
	if (ar.IsShallow())
	{
		if (ar.IsSaving())
		{
			ar << m_mass;
			ar << m_moi;
			ar << m_Fr << m_Mr;
			ar << m_rp << m_rt;
			ar << m_vp << m_vt;
			ar << m_ap << m_at;
			ar << m_qp << m_qt << m_euler;
			ar << m_wp << m_wt;
			ar << m_alp << m_alt;
			for (int i=0; i<6; ++i)
			{
				ar << m_Up[i];
				ar << m_Ut[i];
				ar << m_du[i];
				ar << m_dul[i];
			}
		}
		else
		{
			ar >> m_mass;
			ar >> m_moi;
			ar >> m_Fr >> m_Mr;
			ar >> m_rp >> m_rt;
			ar >> m_vp >> m_vt;
			ar >> m_ap >> m_at;
			ar >> m_qp >> m_qt >> m_euler;
			ar >> m_wp >> m_wt;
			ar >> m_alp >> m_alt;
			for (int i=0; i<6; ++i)
			{
				ar >> m_Up[i];
				ar >> m_Ut[i];
				ar >> m_du[i];
				ar >> m_dul[i];
			}
		}
	}
	else
	{
		ar & m_nID & m_mat & m_mass & m_moi & m_Fr & m_Mr;
		ar & m_Fp & m_Mp;
		ar & m_r0 & m_rp & m_rt & m_vp & m_vt & m_ap & m_at;
		ar & m_qp & m_qt & m_euler & m_wp & m_wt & m_alp & m_alt;
		ar & m_hp & m_ht & m_dhp & m_dht;
		ar & m_bpofr;

		if (ar.IsSaving())
		{
			ar.write(m_BC , sizeof(int), 6);
			ar.write(m_LM , sizeof(int), 6);
			ar.write(m_Up , sizeof(double), 6);
			ar.write(m_Ut , sizeof(double), 6);
			ar.write(m_du , sizeof(double), 6);
			ar.write(m_dul, sizeof(double), 6);
		}
		else
		{
			ar.read(m_BC , sizeof(int), 6);
			ar.read(m_LM , sizeof(int   ), 6);
			ar.read(m_Up , sizeof(double), 6);
			ar.read(m_Ut , sizeof(double), 6);
			ar.read(m_du , sizeof(double), 6);
			ar.read(m_dul, sizeof(double), 6);
		}

		for (int i = 0; i < 6; ++i) ar & m_pDC[i];

		// TODO: should I store parent rigid body (m_prb)
	}
}
