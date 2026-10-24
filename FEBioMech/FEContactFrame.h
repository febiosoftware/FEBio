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

private:
	void Evaluate(FESurface& surf, FESurfaceElement& el, const double* Hr, const double* Hs, bool prev)
	{
		FEMesh& mesh = *surf.GetMesh();
		bool bshb = surf.IsShellBottom();
		vec3d g1(0,0,0), g2(0,0,0);
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
		Set(g1, g2);
	}
};
