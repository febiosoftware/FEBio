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



#pragma once
#include <FECore/FESurfaceLoad.h>
#include <FECore/FEModelParam.h>
#include <vector>

//-----------------------------------------------------------------------------
//! Objective (frame-indifferent) Robin-type external support condition.
//!
//! The surface is attached to a moving anchor (surrounding tissue) by a
//! distribution of anisotropic springs and dashpots. Per unit REFERENCE area:
//!
//!     t0 = -K g - C gdot - p0 * n * (da/dA)
//!
//!     g    = x - y                    (gap to anchor)
//!     gdot = R d/dt(R^T g)            (objective_rate = 1, anchor-frame rate)
//!          = d/dt(g)                  (objective_rate = 0)
//!     K    = kn n~ (x) n~ + kt (I - n~ (x) n~)
//!     C    = cn n~ (x) n~ + ct (I - n~ (x) n~)
//!     n~   = R(t) N                   (reference normal carried by the anchor)
//!
//! Anchor kinematics:
//!     y(X,t) = Xc + c(t) + R(t) (X - Xc) + w(X,t)
//!
//! with R(t) = exp([theta(t)]) built from a rotation vector theta. Under a
//! superposed rigid motion x* = Q x + b applied to BOTH the wall and the anchor,
//! g, gdot, n~ transform as vectors and K, C as second-order tensors, so the
//! traction transforms as t0* = Q t0 (objective). The classical support
//! condition (k u + c udot) is recovered for a fixed anchor (c = theta = w = 0).
//!
//! Time integration:
//!   The rate is computed from the tracked gap variable Z (Z = R^T g if
//!   objective_rate, else Z = g):
//!   - quasi-static analysis: backward Euler, Zdot = (Z - Zn)/dt
//!   - dynamic analysis: Newmark recurrence with the solver's beta, gamma,
//!     and the spring/dashpot force is evaluated at t_{n+alpha_f}:
//!         f_alpha = alpha_f f_{n+1} + (1 - alpha_f) f_n
//!     which matches FEBio's generalized-alpha scheme (alpha_f = 1 for
//!     Euler/Newmark). Because Z is unchanged by any superposed rigid motion
//!     when objective_rate = 1, the discrete rate is exactly objective.
//!
//!   The follower pressure p0 is evaluated at t_{n+1}, like FEPressureLoad.
//!
class FERobinBC : public FESurfaceLoad
{
public:
	FERobinBC(FEModel* pfem);

	void SetSurface(FESurface* ps) override;

	bool Init() override;

	void Activate() override;

	void Update() override;

	void Serialize(DumpStream& ar) override;

public:
	void LoadVector(FEGlobalVector& R) override;

	void StiffnessMatrix(FELinearSystem& LS) override;

public:
	//! anchor rotation at the current time
	mat3d AnchorRotation() const;

	//! anchor position at a surface integration point (current parameter values)
	vec3d AnchorPosition(FESurfaceMaterialPoint& pt, const vec3d& X);

private:
	struct Kinematics
	{
		vec3d	X;		// reference position
		vec3d	x;		// current position
		vec3d	xp;		// position at previous time step
		vec3d	N;		// unit reference normal
		double	J0;		// reference area Jacobian
	};

	void EvalKinematics(FESurfaceElement& el, int n, Kinematics& k);

	mat3ds NormalProjector(const mat3d& Ra, const Kinematics& k) const { return dyad(Ra*k.N); }

	// gap state of the support at one integration point
	struct GapState
	{
		vec3d	Z;		// tracked gap (anchor frame if objective rate)
		vec3d	Zd;		// its rate
		vec3d	Zdd;	// its second rate (Newmark)
		vec3d	f;		// spring/dashpot force per reference area (spatial)
	};

	// make sure the anchor and committed gap state are consistent with the current time
	void SyncAnchor();
	void EvalAnchor(std::vector<vec3d>& y);
	void InitCommitted();

	// evaluate the trial gap state and force at t_{n+1}
	void EvalTrial(FESurfaceMaterialPoint& pt, const Kinematics& k, const mat3d& Ra, GapState& s);

	// d(Zdot)/d(Z) for the active integrator
	double RateFactor() const;

	bool IsDynamic() const;

	int PointIndex(FESurfaceMaterialPoint& pt);

protected:
	FEParamDouble	m_kn;		//!< normal spring stiffness (per reference area)
	FEParamDouble	m_kt;		//!< tangential spring stiffness
	FEParamDouble	m_cn;		//!< normal dashpot coefficient
	FEParamDouble	m_ct;		//!< tangential dashpot coefficient
	FEParamDouble	m_p0;		//!< external (follower) pressure

	vec3d			m_Xc;		//!< anchor center of rotation (reference)
	vec3d			m_c;		//!< anchor translation c(t)
	vec3d			m_theta;	//!< anchor rotation vector theta(t) [rad]
	FEParamVec3		m_w;		//!< additional (mapped) anchor displacement w(X,t)

	bool			m_bshellb;
	bool			m_bobjrate;	//!< use the anchor-frame (objective) rate

private:
	// Committed state (last converged time m_tc) and trial state (last evaluation,
	// at time m_ttr). The trial state is committed when the time advances past
	// it, so failed increments (cutbacks) are discarded automatically.
	std::vector<int>		m_offset;	// first integration point index of each element
	std::vector<vec3d>		m_yt;		// anchor position at m_ty
	std::vector<GapState>	m_com;		// committed state
	std::vector<GapState>	m_trl;		// trial state
	double					m_ty;		// time of m_yt
	double					m_tc;		// time of committed state
	double					m_ttr;		// time of trial state
	bool					m_bcom;		// committed state initialized?
	bool					m_bwarn;

	DECLARE_FECORE_CLASS();
};
