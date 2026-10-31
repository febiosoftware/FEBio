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
#include "FERobinBC.h"
#include "FEBioMech.h"
#include "FESolidAnalysis.h"
#include <FECore/FEFacetSet.h>
#include <FECore/FEModel.h>
#include <FECore/FEAnalysis.h>
#include <FECore/DumpStream.h>
#include <FECore/quatd.h>
#include <FECore/log.h>
#include <cmath>

//=============================================================================
BEGIN_FECORE_CLASS(FERobinBC, FESurfaceLoad)
	ADD_PARAMETER(m_kn      , "normal_stiffness"    )->setUnits("P/L");
	ADD_PARAMETER(m_kt      , "tangential_stiffness")->setUnits("P/L");
	ADD_PARAMETER(m_cn      , "normal_damping"      )->setUnits("P.t/L");
	ADD_PARAMETER(m_ct      , "tangential_damping"  )->setUnits("P.t/L");
	ADD_PARAMETER(m_p0      , "external_pressure"   )->setUnits(UNIT_PRESSURE);
	ADD_PARAMETER(m_Xc      , "anchor_center"       )->setUnits(UNIT_LENGTH);
	ADD_PARAMETER(m_c       , "anchor_translation"  )->setUnits(UNIT_LENGTH);
	ADD_PARAMETER(m_theta   , "anchor_rotation"     );
	ADD_PARAMETER(m_w       , "anchor_displacement" )->setUnits(UNIT_LENGTH);
	ADD_PARAMETER(m_bobjrate, "objective_rate"      );
	ADD_PARAMETER(m_bshellb , "shell_bottom"        );
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FERobinBC::FERobinBC(FEModel* pfem) : FESurfaceLoad(pfem)
{
	m_kn = 0.0;
	m_kt = 0.0;
	m_cn = 0.0;
	m_ct = 0.0;
	m_p0 = 0.0;
	m_Xc    = vec3d(0, 0, 0);
	m_c     = vec3d(0, 0, 0);
	m_theta = vec3d(0, 0, 0);
	m_w     = vec3d(0, 0, 0);
	m_bshellb  = false;
	m_bobjrate = true;

	m_ty = m_tc = m_ttr = -1e99;
	m_bcom  = false;
	m_bwarn = false;
}

//-----------------------------------------------------------------------------
void FERobinBC::SetSurface(FESurface* ps)
{
	FESurfaceLoad::SetSurface(ps);
	FEFacetSet* fs = ps->GetFacetSet();
	m_kn.SetItemList(fs);
	m_kt.SetItemList(fs);
	m_cn.SetItemList(fs);
	m_ct.SetItemList(fs);
	m_p0.SetItemList(fs);
	m_w .SetItemList(fs);
}

//-----------------------------------------------------------------------------
bool FERobinBC::Init()
{
	FESurface& surf = GetSurface();
	surf.SetShellBottom(m_bshellb);

	m_dof.Clear();
	if (m_bshellb == false)
		m_dof.AddVariable(FEBioMech::GetVariableName(FEBioMech::DISPLACEMENT));
	else
		m_dof.AddVariable(FEBioMech::GetVariableName(FEBioMech::SHELL_DISPLACEMENT));
	if (m_dof.IsEmpty()) return false;

	if (FESurfaceLoad::Init() == false) return false;

	// integration point bookkeeping
	int NE = surf.Elements();
	m_offset.assign(NE + 1, 0);
	for (int i = 0; i < NE; ++i)
	{
		FESurfaceElement& el = surf.Element(i);
		if (el.GetLocalID() != i) return false;
		m_offset[i + 1] = m_offset[i] + el.GaussPoints();
	}
	int NP = m_offset[NE];

	GapState s0;
	s0.Z = s0.Zd = s0.Zdd = s0.f = vec3d(0, 0, 0);
	m_yt.assign(NP, vec3d(0, 0, 0));
	m_com.assign(NP, s0);
	m_trl.assign(NP, s0);
	m_ty = m_tc = m_ttr = -1e99;
	m_bcom = false;

	return true;
}

//-----------------------------------------------------------------------------
// Load parameters are scaled by their load controllers only after Init(), so
// the initial anchor configuration and gap state are captured on activation.
void FERobinBC::Activate()
{
	FESurfaceLoad::Activate();
	SyncAnchor();
}

//-----------------------------------------------------------------------------
void FERobinBC::Update()
{
	SyncAnchor();
}

//-----------------------------------------------------------------------------
bool FERobinBC::IsDynamic() const
{
	FEAnalysis* step = GetFEModel()->GetCurrentStep();
	return (step && (step->m_nanalysis == FESolidAnalysis::DYNAMIC));
}

//-----------------------------------------------------------------------------
double FERobinBC::RateFactor() const
{
	const FETimeInfo& tp = GetFEModel()->GetTime();
	double dt = tp.timeIncrement;
	if (dt <= 0) return 0.0;
	if (IsDynamic()) return tp.gamma/(tp.beta*dt);
	return 1.0/dt;
}

//-----------------------------------------------------------------------------
mat3d FERobinBC::AnchorRotation() const
{
	vec3d th = m_theta;
	if (th.norm() < 1e-15) return mat3d::identity();
	quatd q(th);
	return q.RotationMatrix();
}

//-----------------------------------------------------------------------------
vec3d FERobinBC::AnchorPosition(FESurfaceMaterialPoint& pt, const vec3d& X)
{
	mat3d R = AnchorRotation();
	return m_Xc + m_c + R*(X - m_Xc) + m_w(pt);
}

//-----------------------------------------------------------------------------
int FERobinBC::PointIndex(FESurfaceMaterialPoint& pt)
{
	return m_offset[pt.SurfaceElement()->GetLocalID()] + pt.m_index;
}

//-----------------------------------------------------------------------------
void FERobinBC::EvalKinematics(FESurfaceElement& el, int n, Kinematics& k)
{
	FESurface& surf = GetSurface();
	vec3d r0[FEElement::MAX_NODES], rt[FEElement::MAX_NODES], rp[FEElement::MAX_NODES];
	surf.GetReferenceNodalCoordinates(el, r0);
	surf.GetNodalCoordinates(el, rt);
	surf.PreviousNodalCoordinates(el, rp);

	double* H  = el.H(n);
	double* Gr = el.Gr(n);
	double* Gs = el.Gs(n);
	int neln = el.Nodes();

	vec3d G1(0, 0, 0), G2(0, 0, 0);
	k.X = k.x = k.xp = vec3d(0, 0, 0);
	for (int i = 0; i < neln; ++i)
	{
		k.X  += r0[i]*H[i];
		k.x  += rt[i]*H[i];
		k.xp += rp[i]*H[i];
		G1 += r0[i]*Gr[i];
		G2 += r0[i]*Gs[i];
	}
	k.N = G1 ^ G2;
	k.J0 = k.N.unit();
}

//-----------------------------------------------------------------------------
void FERobinBC::EvalAnchor(std::vector<vec3d>& y)
{
	FESurface& surf = GetSurface();
	mat3d R = AnchorRotation();
	y.resize(m_offset.back());
	for (int i = 0; i < surf.Elements(); ++i)
	{
		FESurfaceElement& el = surf.Element(i);
		for (int n = 0; n < el.GaussPoints(); ++n)
		{
			FESurfaceMaterialPoint& pt = static_cast<FESurfaceMaterialPoint&>(*el.GetMaterialPoint(n));
			Kinematics k;
			EvalKinematics(el, n, k);
			y[m_offset[i] + n] = m_Xc + m_c + R*(k.X - m_Xc) + m_w(pt);
		}
	}
}

//-----------------------------------------------------------------------------
// Initialize the committed state from the current configuration. The support
// is assumed to start at rest relative to the anchor (Zdot = Zddot = 0).
void FERobinBC::InitCommitted()
{
	FESurface& surf = GetSurface();
	mat3d Ra = AnchorRotation();
	mat3d RaT = Ra.transpose();
	for (int i = 0; i < surf.Elements(); ++i)
	{
		FESurfaceElement& el = surf.Element(i);
		for (int n = 0; n < el.GaussPoints(); ++n)
		{
			FESurfaceMaterialPoint& pt = static_cast<FESurfaceMaterialPoint&>(*el.GetMaterialPoint(n));
			Kinematics k;
			EvalKinematics(el, n, k);
			int ip = m_offset[i] + n;

			vec3d g = k.x - m_yt[ip];
			mat3ds P = NormalProjector(Ra, k);
			mat3ds Q = mat3ds(1.0, 1.0, 1.0, 0.0, 0.0, 0.0) - P;
			mat3ds K = P*m_kn(pt) + Q*m_kt(pt);

			GapState& s = m_com[ip];
			s.Z   = (m_bobjrate ? RaT*g : g);
			s.Zd  = vec3d(0, 0, 0);
			s.Zdd = vec3d(0, 0, 0);
			s.f   = -(K*g);
		}
	}
	m_trl = m_com;
	m_tc = m_ttr = GetFEModel()->GetTime().currentTime;
	m_bcom = true;
}

//-----------------------------------------------------------------------------
void FERobinBC::SyncAnchor()
{
	if (m_offset.empty()) return;

	const FETimeInfo& tp = GetFEModel()->GetTime();
	double t  = tp.currentTime;
	double dt = tp.timeIncrement;
	double tol = 1e-12*(1.0 + fabs(t)) + 1e-6*fabs(dt);

	// anchor at the current time (anchor parameters depend on time only)
	if (fabs(t - m_ty) > tol)
	{
		EvalAnchor(m_yt);
		m_ty = t;
	}

	// first call: initialize the committed state
	if (!m_bcom) { InitCommitted(); return; }

	// Time has advanced past the committed state: the last trial state was the
	// converged solution of the previous increment, so commit it. After a
	// cutback the previous time equals m_tc and nothing is committed.
	if (dt > 0)
	{
		double tprev = t - dt;
		if (fabs(tprev - m_tc) > tol)
		{
			if (fabs(m_ttr - tprev) <= tol)
			{
				m_com = m_trl;
				m_tc = m_ttr;
			}
			else
			{
				if (!m_bwarn)
				{
					feLogWarning("Robin BC: no converged state found at t = %lg; keeping the last committed state.", tprev);
					m_bwarn = true;
				}
				m_tc = tprev;
			}
		}
	}
}

//-----------------------------------------------------------------------------
// Trial gap state and spring/dashpot force (per reference area) at t_{n+1}.
void FERobinBC::EvalTrial(FESurfaceMaterialPoint& pt, const Kinematics& k, const mat3d& Ra, GapState& s)
{
	const FETimeInfo& tp = GetFEModel()->GetTime();
	double dt = tp.timeIncrement;

	int ip = PointIndex(pt);
	const GapState& c = m_com[ip];

	// gap and tracked variable (anchor frame for the objective rate)
	vec3d g = k.x - m_yt[ip];
	s.Z = (m_bobjrate ? Ra.transpose()*g : g);

	// rate of the tracked variable
	if (dt <= 0)
	{
		s.Zd = s.Zdd = vec3d(0, 0, 0);
	}
	else if (IsDynamic())
	{
		double beta = tp.beta, gamma = tp.gamma;
		s.Zdd = (s.Z - c.Z)/(beta*dt*dt) - c.Zd/(beta*dt) - c.Zdd*(0.5/beta - 1.0);
		s.Zd  = c.Zd + (c.Zdd*(1.0 - gamma) + s.Zdd*gamma)*dt;
	}
	else
	{
		s.Zd  = (s.Z - c.Z)/dt;
		s.Zdd = vec3d(0, 0, 0);
	}

	// spatial gap rate
	vec3d gdot = (m_bobjrate ? Ra*s.Zd : s.Zd);

	// anchor-convected normal and support tensors
	mat3ds P = NormalProjector(Ra, k);
	mat3ds Q = mat3ds(1.0, 1.0, 1.0, 0.0, 0.0, 0.0) - P;
	mat3ds K = P*m_kn(pt) + Q*m_kt(pt);
	mat3ds C = P*m_cn(pt) + Q*m_ct(pt);

	s.f = -(K*g + C*gdot);
}

//-----------------------------------------------------------------------------
void FERobinBC::LoadVector(FEGlobalVector& R)
{
	const FETimeInfo& tp = GetFEModel()->GetTime();
	if (tp.currentTime == 0) return;

	SyncAnchor();
	mat3d Ra = AnchorRotation();
	double alpha = (IsDynamic() ? tp.alphaf : 1.0);

	FESurface& surf = GetSurface();
	surf.LoadVector(R, m_dof, false, [&](FESurfaceMaterialPoint& pt, const FESurfaceDofShape& dof_a, std::vector<double>& val) {

		FESurfaceElement& el = *pt.SurfaceElement();
		Kinematics k;
		EvalKinematics(el, pt.m_index, k);

		// trial state at t_{n+1} (each integration point belongs to one
		// element, so this write is thread-safe)
		int ip = PointIndex(pt);
		GapState s;
		EvalTrial(pt, k, Ra, s);
		m_trl[ip] = s;

		// spring/dashpot force at t_{n+alpha_f}
		vec3d f = (s.f*alpha + m_com[ip].f*(1.0 - alpha))*k.J0;

		// follower external pressure on current area: -p0 n da
		double p0 = m_p0(pt);
		f -= (pt.dxr ^ pt.dxs)*p0;

		if (m_bshellb) f = -f;

		double Na = dof_a.shape;
		val[0] = Na*f.x;
		val[1] = Na*f.y;
		val[2] = Na*f.z;
	});
	m_ttr = tp.currentTime;
}

//-----------------------------------------------------------------------------
void FERobinBC::StiffnessMatrix(FELinearSystem& LS)
{
	const FETimeInfo& tp = GetFEModel()->GetTime();
	if (tp.currentTime == 0) return;

	mat3d Ra = AnchorRotation();
	double alpha = (IsDynamic() ? tp.alphaf : 1.0);
	double r = RateFactor();	// d(gdot)/d(x)

	FESurface& surf = GetSurface();
	surf.LoadStiffness(LS, m_dof, m_dof, [&](FESurfaceMaterialPoint& pt, const FESurfaceDofShape& dof_a, const FESurfaceDofShape& dof_b, matrix& kab) {

		FESurfaceElement& el = *pt.SurfaceElement();
		Kinematics k;
		EvalKinematics(el, pt.m_index, k);

		mat3ds P = NormalProjector(Ra, k);
		mat3ds Q = mat3ds(1.0, 1.0, 1.0, 0.0, 0.0, 0.0) - P;

		double kn = m_kn(pt), kt = m_kt(pt);
		double cn = m_cn(pt), ct = m_ct(pt);

		// d(f_alpha)/d(x_{n+1}) = -alpha (K + r C); with the objective rate
		// d(gdot)/dx = R r R^T = r I, so the form is the same in both cases
		mat3ds Keff = (P*(kn + cn*r) + Q*(kt + ct*r))*alpha;

		double Na = dof_a.shape;
		double Nb = dof_b.shape;

		// spring/dashpot part (reference area: symmetric, no area variation)
		mat3d Kab = Keff*(Na*Nb*k.J0);

		// follower pressure part (same form as FEPressureLoad, non-symmetric)
		double P0 = -m_p0(pt);
		if (m_bshellb) P0 = -P0;
		vec3d vab = (pt.dxs*dof_b.shape_deriv_r - pt.dxr*dof_b.shape_deriv_s)*(P0*Na);
		Kab += mat3d(mat3da(vab));

		kab.set(0, 0, Kab);
	});
}

//-----------------------------------------------------------------------------
static void SerializeStates(DumpStream& ar, std::vector<vec3d>& Z, std::vector<vec3d>& Zd, std::vector<vec3d>& Zdd, std::vector<vec3d>& f)
{
	ar & Z & Zd & Zdd & f;
}

void FERobinBC::Serialize(DumpStream& ar)
{
	FESurfaceLoad::Serialize(ar);
	if (ar.IsShallow()) return;

	ar & m_offset & m_yt & m_ty & m_tc & m_ttr & m_bcom & m_bwarn;

	std::vector<GapState>* states[2] = { &m_com, &m_trl };
	for (std::vector<GapState>* ps : states)
	{
		std::vector<GapState>& S = *ps;
		std::vector<vec3d> Z, Zd, Zdd, f;
		if (ar.IsSaving())
		{
			for (GapState& s : S) { Z.push_back(s.Z); Zd.push_back(s.Zd); Zdd.push_back(s.Zdd); f.push_back(s.f); }
			SerializeStates(ar, Z, Zd, Zdd, f);
		}
		else
		{
			SerializeStates(ar, Z, Zd, Zdd, f);
			S.resize(Z.size());
			for (size_t i = 0; i < S.size(); ++i) { S[i].Z = Z[i]; S[i].Zd = Zd[i]; S[i].Zdd = Zdd[i]; S[i].f = f[i]; }
		}
	}
}
