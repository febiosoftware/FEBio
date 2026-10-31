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
#include "mat3d.h"
#include "quatd.h"
#include <math.h>

//-----------------------------------------------------------------------------
// Incrementally objective kinematics for the generalized-alpha method.
//
// The generalized-alpha method evaluates the balance equations at the
// intermediate time t_{n+af}. Interpolating the deformation gradient linearly,
// F_{n+af} = af*F_{n+1} + (1-af)*F_n, is not incrementally objective: for a
// rigid-body rotation increment, F_{n+af} is not a rotation, which produces a
// spurious strain of order af*(1-af)*dtheta^2 and a spurious contraction of the
// intermediate configuration. Similarly, the velocity gradient
// L = (F_{n+1} - F_n)/dt . F_{n+af}^{-1} has a non-zero symmetric part (a
// spurious rate of deformation, proportional to (af - 1/2)*dtheta/dt) under
// rigid-body rotation, unless af = 1/2.
//
// This file provides objective alternatives:
// - FEObjectiveInterpolation interpolates the right stretch tensor linearly and
//   the rotation along the geodesic (constant angular velocity), so that a
//   rigid-body increment produces no strain at the intermediate time. It also
//   provides the derivative of the interpolated deformation gradient with
//   respect to F_{n+1}, needed for a consistent tangent.
// - FEMidpointVelocityGradient evaluates L on the midpoint configuration
//   (Hughes and Winget, 1980), whose symmetric part vanishes exactly for a
//   rigid-body increment of any size; this is also the second-order accurate
//   approximation of the velocity gradient over the time step.
//-----------------------------------------------------------------------------

//! axial vector of the skew part of A, i.e. of (A - A^T), such that (A - A^T).b = axial x b
inline vec3d FEAxialVector(const mat3d& A)
{
	return vec3d(A(2,1) - A(1,2), A(0,2) - A(2,0), A(1,0) - A(0,1));
}

//! skew tensor [a]x such that [a]x.b = a x b
inline mat3d FESkew(const vec3d& a)
{
	return mat3d(0.0, -a.z, a.y,
		         a.z, 0.0, -a.x,
		        -a.y, a.x, 0.0);
}

//! rotation tensor from a rotation vector (Rodrigues formula)
inline mat3d FERotationFromVector(const vec3d& phi)
{
	double th = phi.norm();
	if (th < 1e-14) return mat3d(mat3dd(1.0)) + FESkew(phi);
	mat3d K = FESkew(phi / th);
	return mat3d(mat3dd(1.0)) + K*sin(th) + (K*K)*(1.0 - cos(th));
}

//! rotation vector of a rotation tensor (principal value, angle in [0, pi])
inline vec3d FERotationVector(const mat3d& R)
{
	quatd q(R);
	double w = q.w;
	vec3d v(q.x, q.y, q.z);
	if (w < 0) { w = -w; v = -v; }
	double s = v.norm();
	if (s < 1e-14) return v*2.0;
	double th = 2.0*atan2(s, w);
	return v*(th / s);
}

//! right Jacobian of SO(3): exp(phi)^T . d(exp(phi)) = [Jr(phi).dphi]x
inline mat3d FERightJacobianSO3(const vec3d& phi)
{
	double th = phi.norm();
	mat3d K = FESkew(phi);
	mat3d I(mat3dd(1.0));
	if (th < 1e-6) return I - K*0.5 + (K*K)*(1.0/6.0);
	return I - K*((1.0 - cos(th))/(th*th)) + (K*K)*((th - sin(th))/(th*th*th));
}

//! inverse of the right Jacobian of SO(3)
inline mat3d FEInvRightJacobianSO3(const vec3d& phi)
{
	double th = phi.norm();
	mat3d K = FESkew(phi);
	mat3d I(mat3dd(1.0));
	if (th < 1e-6) return I + K*0.5 + (K*K)*(1.0/12.0);
	return I + K*0.5 + (K*K)*(1.0/(th*th) - (1.0 + cos(th))/(2.0*th*sin(th)));
}

//-----------------------------------------------------------------------------
//! Objective interpolation of the deformation gradient at t_{n+a}:
//!     F_a = R_a.U_a,   U_a = a*U_{n+1} + (1-a)*U_n,   R_a = R_n.exp(a*log(R_n^T.R_{n+1})),
//! where F = R.U is the right polar decomposition. For a = 1, F_a = F_{n+1}.
class FEObjectiveInterpolation
{
public:
	FEObjectiveInterpolation() : m_a(1.0) {}

	//! initialize from the deformation gradients at t_n (Fp) and t_{n+1} (Ft)
	void Init(const mat3d& Fp, const mat3d& Ft, double a)
	{
		m_a = a;
		if (a == 1.0) { m_Fa = Ft; return; }

		mat3d Rp; mat3ds Up;
		Fp.right_polar(Rp, Up);
		Ft.right_polar(m_Rt, m_Ut);
		vec3d phi = FERotationVector(Rp.transpose()*m_Rt);
		m_Ra = Rp*FERotationFromVector(phi*a);
		m_Ua = m_Ut*a + Up*(1.0 - a);
		m_Fa = m_Ra*mat3d(m_Ua);

		// quantities needed for the derivative
		m_Mi = (mat3d(mat3dd(m_Ut.tr())) - mat3d(m_Ut)).inverse();
		m_T = FERightJacobianSO3(phi*a)*FEInvRightJacobianSO3(phi)*a;
	}

	//! the interpolated deformation gradient
	const mat3d& F() const { return m_Fa; }

	//! interpolation factor
	double alpha() const { return m_a; }

	//! derivative of F_a in the direction dFt of F_{n+1}
	mat3d dF(const mat3d& dFt) const
	{
		if (m_a == 1.0) return dFt;

		// derivative of the polar decomposition of F_{n+1}:
		// dR = R.[w]x, dU = R^T.dF - [w]x.U, (tr(U) I - U).w = axial(R^T.dF - dF^T.R)
		mat3d A = m_Rt.transpose()*dFt;
		vec3d w = m_Mi*FEAxialVector(A);
		mat3d dU = A - FESkew(w)*mat3d(m_Ut);

		// derivative of R_a = R_n.exp(a*phi): R_a^T.dR_a = [a Jr(a phi).Jr(phi)^-1.w]x
		vec3d wa = m_T*w;
		return m_Ra*(FESkew(wa)*mat3d(m_Ua) + dU*m_a);
	}

private:
	double	m_a;
	mat3d	m_Fa;
	mat3d	m_Rt, m_Ra;
	mat3ds	m_Ut, m_Ua;
	mat3d	m_Mi;	// (tr(U_{n+1}) I - U_{n+1})^-1
	mat3d	m_T;	// a*Jr(a phi).Jr(phi)^-1
};

//! Objective interpolation of the deformation gradient (see FEObjectiveInterpolation)
inline mat3d FEObjectiveInterpolateF(const mat3d& Fp, const mat3d& Ft, double a)
{
	FEObjectiveInterpolation oi;
	oi.Init(Fp, Ft, a);
	return oi.F();
}

//! Velocity gradient over the time step, evaluated on the midpoint configuration:
//! L = (F_{n+1} - F_n)/dt . F_{n+1/2}^{-1}, F_{n+1/2} = (F_{n+1} + F_n)/2.
inline mat3d FEMidpointVelocityGradient(const mat3d& Fp, const mat3d& Ft, double dt)
{
	mat3d Fm = (Ft + Fp)*0.5;
	return (Ft - Fp)*Fm.inverse()*(1.0 / dt);
}

//! Derivative of the midpoint velocity gradient in the direction dFt of F_{n+1}:
//! dL = (I/dt - L/2).dFt.F_{n+1/2}^{-1}
inline mat3d FEMidpointVelocityGradientDerivative(const mat3d& Fp, const mat3d& Ft, double dt, const mat3d& dFt)
{
	mat3d Fmi = ((Ft + Fp)*0.5).inverse();
	mat3d L = (Ft - Fp)*Fmi*(1.0 / dt);
	return (mat3d(mat3dd(1.0/dt)) - L*0.5)*dFt*Fmi;
}
