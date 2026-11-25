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
#include <FECore/FESurface.h>
#include <FECore/FEMesh.h>
#include <FECore/DumpStream.h>
#include <FECore/mat3d.h>

//-----------------------------------------------------------------------------
//! Orthonormal corotational frame R = [e1 e2 n] attached to a material point of
//! a contact surface: e1 = g1/|g1|, n = g1 x g2/|g1 x g2|, e2 = n x e1, where
//! g1, g2 are the covariant base vectors of a surface element.
//!
//! Contact interfaces that store spatial vectors across time steps (vectorial
//! Lagrange multipliers, initial gap vectors) use this frame to push those
//! vectors forward with the local rotation of the interface at the start of
//! each time step (see FEContactFrame::PushForward). Without this, the stored
//! vectors remain fixed in space while the interface rotates, which violates
//! frame invariance (objectivity).
//!
//! NOTE: Pushing forward at the start of each time step is adequate for quantities
//! that are re-established by the solution itself (e.g. Lagrange multipliers, which
//! augmentation drives to the current traction). It is NOT adequate for a stored
//! geometric offset such as an initial gap vector Gap: holding R_n.Gap fixed within
//! a step while the interface rotates by dR produces a spurious gap (dR - I).R_n.Gap,
//! i.e. spurious tractions ~ eps*|dtheta_step|*|Gap| that depend on the time step
//! size. Such offsets must be stored by their components in the frame and rotated
//! with the frame evaluated in the CURRENT configuration (prev = false). Since R
//! then depends on the nodal positions, the tie also exerts forces through R, which
//! are needed for moment balance (see FETiedFrameOffset at the end of this file).
class FEContactFrame
{
public:
	vec3d	e1, e2, n;

public:
	FEContactFrame() : e1(0,0,0), e2(0,0,0), n(0,0,0) {}

	//! true if the frame has been initialized
	bool IsValid() const { return (n*n > 0.0); }

	//! R.a
	vec3d R(const vec3d& a) const { return e1*a.x + e2*a.y + n*a.z; }

	//! R^T.a
	vec3d Rt(const vec3d& a) const { return vec3d(e1*a, e2*a, n*a); }

	//! build the frame from two covariant base vectors
	void Set(const vec3d& g1, const vec3d& g2)
	{
		n = g1 ^ g2; n.unit();
		e1 = g1; e1.unit();
		e2 = n ^ e1;
	}

	//! Evaluate the frame at natural coordinates (r,s) of a surface element.
	//! If prev is true, the configuration at the start of the time step (m_rp) is used,
	//! otherwise the current configuration (m_rt).
	void Evaluate(FESurface& surf, FESurfaceElement& el, double r, double s, bool prev = true)
	{
		double Hr[FEElement::MAX_NODES], Hs[FEElement::MAX_NODES];
		el.shape_deriv(Hr, Hs, r, s);
		Evaluate(surf, el, Hr, Hs, prev);
	}

	//! Evaluate the frame at integration point j of a surface element
	void Evaluate(FESurface& surf, FESurfaceElement& el, int j, bool prev = true)
	{
		Evaluate(surf, el, el.Gr(j), el.Gs(j), prev);
	}

	//! Evaluate the frame at the center of a surface element
	void EvaluateCenter(FESurface& surf, FESurfaceElement& el, bool prev = true)
	{
		int ne = el.Nodes();
		bool btri = ((ne == 3) || (ne == 6) || (ne == 7));
		double c = (btri ? 1.0/3.0 : 0.0);
		Evaluate(surf, el, c, c, prev);
	}

	//! Rotate a vector v, defined relative to the frame Rold, so that it is
	//! defined relative to the frame Rnew: v <- Rnew.Rold^T.v
	static vec3d PushForward(const FEContactFrame& Rold, const FEContactFrame& Rnew, const vec3d& v)
	{
		return Rnew.R(Rold.Rt(v));
	}

	//! Update a stored frame to Rnew and push the vectors v[0..nv-1] forward
	//! with the incremental rotation. If the stored frame was not yet
	//! initialized, it is simply set to Rnew.
	void Update(const FEContactFrame& Rnew, vec3d* v, int nv)
	{
		if (IsValid())
		{
			for (int i = 0; i < nv; ++i) v[i] = PushForward(*this, Rnew, v[i]);
		}
		*this = Rnew;
	}

	void Serialize(DumpStream& ar)
	{
		ar & e1 & e2 & n;
	}

	//! Evaluate the frame at integration point j of a surface element, and also
	//! return the generating vectors g1, g2 that were used (including the sign
	//! flip for shell-bottom surfaces), as needed by Derivatives().
	//! The derivative of g1 (g2) with respect to the position of node c is
	//! sgn*Gr[c] (sgn*Gs[c]), where sgn = -1 for shell-bottom surfaces, +1 otherwise.
	void Evaluate(FESurface& surf, FESurfaceElement& el, int j, bool prev, vec3d& g1, vec3d& g2, double& sgn)
	{
		BaseVectors(surf, el, el.Gr(j), el.Gs(j), prev, g1, g2);
		sgn = (surf.IsShellBottom() ? -1.0 : 1.0);
		Set(g1, g2);
	}

	//! Linearization of the rotated vector R.v, for fixed frame components v, with
	//! respect to the generating vectors of the frame:
	//!     d(R.v) = Da.dg1 + Db.dg2
	//! The frame must have been built (Set) from the same g1 and g2.
	void Derivatives(const vec3d& g1, const vec3d& g2, const vec3d& v, mat3d& Da, mat3d& Db) const
	{
		mat3d I; I.unit();
		double la = g1.norm();
		double lm = (g1 ^ g2).norm();

		mat3d P1 = (I - (e1 & e1))/la;			// d(e1)/d(g1)
		mat3d Pn = (I - (n & n))/lm;			// d(n)/d(g1 x g2)
		mat3d dn_da = -(Pn*Skew(g2));			// d(g1 x g2)/d(g1) = -[g2]x
		mat3d dn_db =   Pn*Skew(g1);			// d(g1 x g2)/d(g2) =  [g1]x
		mat3d E1 = Skew(e1), N = Skew(n);
		mat3d de2_da = -(E1*dn_da) + N*P1;		// e2 = n x e1
		mat3d de2_db = -(E1*dn_db);

		Da = P1*v.x + de2_da*v.y + dn_da*v.z;
		Db =          de2_db*v.y + dn_db*v.z;
	}

	//! skew-symmetric matrix [w]x such that [w]x.a = w x a
	static mat3d Skew(const vec3d& w)
	{
		return mat3d(0.0, -w.z,  w.y,
		             w.z,  0.0, -w.x,
		            -w.y,  w.x,  0.0);
	}

	//! Second derivatives of the scalar s.R(g1,g2).v (s and v fixed) with respect to
	//! the generating vectors g = (g1, g2), returned as 3x3 blocks
	//!     H[k][l] = d2(s.R.v)/(dg_k dg_l),   k,l = 0 (g1), 1 (g2),
	//! so that H[1][0] = H[0][1]^T. The frame must have been Set from g1 and g2.
	//! With R = [e1 e2 n], s.R.v = v.x (s.e1) + v.y (s.e2) + v.z (s.n), where
	//! e1 = g1/|g1|, n = m/|m| with m = g1 x g2, and s.e2 = s.(n x e1) = n.(e1 x s).
	void Hessian(const vec3d& g1, const vec3d& g2, const vec3d& v, const vec3d& s, mat3d H[2][2]) const
	{
		mat3d I; I.unit();
		double la = g1.norm();
		double lm = (g1 ^ g2).norm();
		mat3d P1 = I - (e1 & e1);
		mat3d Pn = I - (n & n);
		mat3d A[2] = { -Skew(g2), Skew(g1) };			// dm/dg1, dm/dg2
		mat3d Je1 = P1/la;								// de1/dg1 (de1/dg2 = 0)
		mat3d Jn[2] = { (Pn/lm)*A[0], (Pn/lm)*A[1] };	// dn/dg1, dn/dg2

		for (int k = 0; k < 2; ++k) for (int l = 0; l < 2; ++l) H[k][l].zero();

		// Hessian of a.e1 (a fixed): only the g1-g1 block
		auto addE1 = [&](const vec3d& a, double c) {
			if (c == 0.0) return;
			H[0][0] += UnitHessian(la, e1, P1, a)*c;
		};
		// Hessian of a.n (a fixed), including the second derivative of m = g1 x g2
		auto addN = [&](const vec3d& a, double c) {
			if (c == 0.0) return;
			mat3d Hm = UnitHessian(lm, n, Pn, a);
			vec3d q = (Pn*a)/lm;						// d(a.n)/dm
			mat3d H01 = A[0].transpose()*Hm*A[1] - Skew(q);
			H[0][0] += (A[0].transpose()*Hm*A[0])*c;
			H[1][1] += (A[1].transpose()*Hm*A[1])*c;
			H[0][1] += H01*c;
			H[1][0] += H01.transpose()*c;
		};

		// v.x (s.e1)
		addE1(s, v.x);

		// v.z (s.n)
		addN(s, v.z);

		// v.y (s.e2), with s.e2 = n^T.B.e1, B = -[s]x
		if (v.y != 0.0)
		{
			addN(e1 ^ s, v.y);							// d(s.e2)/dn  = e1 x s
			addE1(s ^ n, v.y);							// d(s.e2)/de1 = s x n
			mat3d B = -Skew(s);
			mat3d C0 = Jn[0].transpose()*B*Je1;			// cross terms Jn^T.B.Je1 + Je1^T.B^T.Jn
			mat3d C1 = Jn[1].transpose()*B*Je1;
			H[0][0] += (C0 + C0.transpose())*v.y;
			H[1][0] += C1*v.y;
			H[0][1] += C1.transpose()*v.y;
		}
	}

	//! Assemble the stiffness of a tied integration point whose contact potential
	//! psi(y, g) depends on the gap vector y and on the frame generating vectors
	//! g = (g1, g2) of the primary surface:
	//!     ke(a,b) += scale*[ Y_a Y_b eps I (if bYY) + Y_a sum_k T_b[k] Hyg[k]
	//!                        + sum_k T_a[k] Hyg[k]^T Y_b + sum_kl T_a[k] T_b[l] Hgg[k][l] ]
	//! where d(y)/d(x_a) = Y[a]*I, d(g_k)/d(x_a) = T[a][k]*I (zero for secondary
	//! nodes), Hyg[k] = d2(psi)/(dy dg_k) and Hgg[k][l] = d2(psi)/(dg_k dg_l).
	//! off[a] is the offset of the displacement dofs of node a in ke.
	template <class MAT>
	static void AddHessian(MAT& ke, int nn, const int* off, const double* Y, const double (*T)[2],
		double eps, bool bYY, const mat3d Hyg[2], const mat3d Hgg[2][2], double scale)
	{
		for (int a = 0; a < nn; ++a)
			for (int b = 0; b < nn; ++b)
			{
				mat3d K; K.zero();
				if (bYY) { K(0,0) = K(1,1) = K(2,2) = eps*Y[a]*Y[b]; }
				for (int k = 0; k < 2; ++k)
				{
					if (T[b][k] != 0.0) K += Hyg[k]*(Y[a]*T[b][k]);
					if (T[a][k] != 0.0) K += Hyg[k].transpose()*(T[a][k]*Y[b]);
					for (int l = 0; l < 2; ++l)
						if (T[a][k]*T[b][l] != 0.0) K += Hgg[k][l]*(T[a][k]*T[b][l]);
				}
				for (int l = 0; l < 3; ++l)
					for (int m = 0; m < 3; ++m)
						ke[off[a] + l][off[b] + m] += scale*K(l, m);
			}
	}

private:
	//! Hessian of a.u with respect to x, where u = x/|x| (a fixed), given |x|, u and P = I - u x u
	static mat3d UnitHessian(double lx, const vec3d& u, const mat3d& P, const vec3d& a)
	{
		vec3d Pa = P*a;
		return -(P*(u*a) + (Pa & u) + (u & Pa))/(lx*lx);
	}

	void Evaluate(FESurface& surf, FESurfaceElement& el, const double* Hr, const double* Hs, bool prev)
	{
		vec3d g1, g2;
		BaseVectors(surf, el, Hr, Hs, prev, g1, g2);
		Set(g1, g2);
	}

	static void BaseVectors(FESurface& surf, FESurfaceElement& el, const double* Hr, const double* Hs, bool prev, vec3d& g1, vec3d& g2)
	{
		FEMesh& mesh = *surf.GetMesh();
		bool bshb = surf.IsShellBottom();
		g1 = vec3d(0,0,0); g2 = vec3d(0,0,0);
		int ne = el.Nodes();
		for (int i = 0; i < ne; ++i)
		{
			FENode& nd = mesh.Node(el.m_node[i]);
			vec3d y;
			if (bshb) y = (prev ? nd.sp() : nd.st());
			else      y = (prev ? nd.m_rp : nd.m_rt);
			g1 += y*Hr[i];
			g2 += y*Hs[i];
		}
		if (bshb) { g1 = -g1; g2 = -g2; }
	}
};

//-----------------------------------------------------------------------------
//! Helper for a tied integration point whose initial gap G is stored by its
//! components in the corotational frame R(g1, g2) of the primary surface, so that
//! the gap function is dg = y - R.G, with y the current gap vector (y = x2 - x1,
//! or x1 - x2, depending on the interface's convention).
//!
//! The tie forces derive from the potential (per unit area)
//!     psi(y, g) = lam.dg + eps/2 dg.dg,
//! which, through its dependence on the frame generating vectors g = (g1, g2),
//! also contributes forces on the primary nodes. These forces balance the moment
//! of the traction pair (+t at x1, -t at x2 = x1 + R.G), so that angular momentum
//! is conserved, and they make the formulation variational, which is what keeps
//! the augmented Lagrangian iterations convergent. Omitting them leaves a net
//! couple (R.G) x t at every tied point.
//!
//! Two forms of the multiplier are supported:
//! - EvalSpatial: lam is a spatial vector, held fixed during the iterations;
//! - EvalCorotational: lam = R.L, with L stored by its frame components.
struct FETiedFrameOffset
{
	FEContactFrame	R;			//!< the frame
	mat3d	DG[2];				//!< d(R.G)/dg_k
	vec3d	dpsi_dg[2];			//!< d(psi)/dg_k at fixed y
	mat3d	Hyg[2];				//!< d2(psi)/(dy dg_k)
	mat3d	Hgg[2][2];			//!< d2(psi)/(dg_k dg_l)

	//! lam fixed in space; t = lam + eps*dg is the traction
	void EvalSpatial(const vec3d& g1, const vec3d& g2, const vec3d& G, const vec3d& t, double eps)
	{
		R.Set(g1, g2);
		R.Derivatives(g1, g2, G, DG[0], DG[1]);
		mat3d HG[2][2];
		R.Hessian(g1, g2, G, t, HG);
		for (int k = 0; k < 2; ++k)
		{
			dpsi_dg[k] = -(DG[k].transpose()*t);
			Hyg[k] = -(DG[k]*eps);
			for (int l = 0; l < 2; ++l) Hgg[k][l] = (DG[k].transpose()*DG[l])*eps - HG[k][l];
		}
	}

	//! lam = R.L (L in frame components); y is the gap vector and dg = y - R.G
	void EvalCorotational(const vec3d& g1, const vec3d& g2, const vec3d& G, const vec3d& L,
		const vec3d& y, const vec3d& dg, double eps)
	{
		R.Set(g1, g2);
		R.Derivatives(g1, g2, G, DG[0], DG[1]);
		mat3d DL[2];
		R.Derivatives(g1, g2, L, DL[0], DL[1]);
		mat3d HL[2][2], HG[2][2];
		R.Hessian(g1, g2, L, y, HL);
		R.Hessian(g1, g2, G, dg, HG);
		// psi = (R.L).y + eps/2 |y - R.G|^2 (+ the constant L.G)
		for (int k = 0; k < 2; ++k)
		{
			dpsi_dg[k] = DL[k].transpose()*y - DG[k].transpose()*dg*eps;
			Hyg[k] = DL[k] - DG[k]*eps;
			for (int l = 0; l < 2; ++l) Hgg[k][l] = HL[k][l] - HG[k][l]*eps + (DG[k].transpose()*DG[l])*eps;
		}
	}

	//! frame force on a primary node: -W*sum_k T[k]*d(psi)/dg_k
	vec3d Force(const double T[2], double W) const
	{
		return -(dpsi_dg[0]*T[0] + dpsi_dg[1]*T[1])*W;
	}

	//! For interfaces that integrate over the current area W = J*w: tangent of the
	//! frame forces due to the variation of J,
	//!     ke(a,b) += scale*(sum_k T_a[k] dpsi/dg_k) x dJ/dx_b * w
	//! with dJ/dx_b = (g2 x n) Gr_b + (n x g1) Gs_b (primary nodes, nn1 of them).
	//! This term is non-symmetric.
	template <class MAT>
	void AddAreaTangent(MAT& ke, int nn1, const int* off, const double (*T)[2],
		const vec3d& g1, const vec3d& g2, const double* Gr, const double* Gs, double w, double scale) const
	{
		vec3d nu = g1 ^ g2; nu.unit();
		vec3d c1 = g2 ^ nu, c2 = nu ^ g1;
		for (int a = 0; a < nn1; ++a)
		{
			vec3d fa = dpsi_dg[0]*T[a][0] + dpsi_dg[1]*T[a][1];
			for (int b = 0; b < nn1; ++b)
			{
				vec3d dJ = (c1*Gr[b] + c2*Gs[b])*w;
				for (int l = 0; l < 3; ++l)
					for (int m = 0; m < 3; ++m)
						ke[off[a] + l][off[b] + m] += scale*(&fa.x)[l]*(&dJ.x)[m];
			}
		}
	}
};
