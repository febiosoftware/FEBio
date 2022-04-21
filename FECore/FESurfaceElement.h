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
#include "FEElement.h"

//-----------------------------------------------------------------------------
//!  This class defines a surface element

class FECORE_API FESurfaceElement : public FEElement
{
public:
	FESurfaceElement();

	FESurfaceElement(const FESurfaceElement& el);

	FESurfaceElement& operator = (const FESurfaceElement& el);

	virtual void SetTraits(FEElementTraits* pt) override;

	double* GaussWeights() { return &((FESurfaceElementTraits*)(m_pT))->gw[0]; }			// weights of integration points
	const double* GaussWeights() const { return &((FESurfaceElementTraits*)(m_pT))->gw[0]; }			// weights of integration points
	double gr(int n) const { return ((FESurfaceElementTraits*)(m_pT))->gr[n]; }	// integration point coordinate r
	double gs(int n) const { return ((FESurfaceElementTraits*)(m_pT))->gs[n]; }	// integration point coordinate  s
    double cr() const { return ((FESurfaceElementTraits*)(m_pT))->cr; }    // centroid point coordinate r
    double cs() const { return ((FESurfaceElementTraits*)(m_pT))->cs; }    // centroid point coordinate s

	double* Gr(int n) const { return ((FESurfaceElementTraits*)(m_pT))->Gr[n]; }	// shape function derivative to r
	double* Gs(int n) const { return ((FESurfaceElementTraits*)(m_pT))->Gs[n]; }	// shape function derivative to s

	double eval(double* d, int n)
	{
		double* N = H(n);
		int ne = Nodes();
		double a = 0;
		for (int i=0; i<ne; ++i) a += N[i]*d[i];
		return a;
	}

	double eval(int order, double* d, int n)
	{
		double* N = H(order, n);
		int ne = ShapeFunctions(order);
		double a = 0;
		for (int i = 0; i<ne; ++i) a += N[i] * d[i];
		return a;
	}

	double eval(double* d, double r, double s)
	{
		int n = Nodes();
		double H[FEElement::MAX_NODES];
		shape_fnc(H, r, s);
		double a = 0;
		for (int i=0; i<n; ++i) a += H[i]*d[i];
		return a;
	}

	double eval(int order, double* d, double r, double s)
	{
		int n = ShapeFunctions(order);
		double H[FEElement::MAX_NODES];
		shape_fnc(order, H, r, s);
		double a = 0;
		for (int i = 0; i<n; ++i) a += H[i] * d[i];
		return a;
	}

	vec3d eval(vec3d* d, double r, double s)
	{
		int n = Nodes();
		double H[FEElement::MAX_NODES];
		shape_fnc(H, r, s);
		vec3d a(0,0,0);
		for (int i=0; i<n; ++i) a += d[i]*H[i];
		return a;
	}

	vec3d eval(vec3d* d, int n)
	{
		int ne = Nodes();
		double* N = H(n);
		vec3d a(0,0,0);
		for (int i=0; i<ne; ++i) a += d[i]*N[i];
		return a;
	}

	double eval_deriv1(double* d, int j)
	{
		double* Hr = Gr(j);
		int n = Nodes();
		double s = 0;
		for (int i=0; i<n; ++i) s +=  Hr[i]*d[i];
		return s;
	}

	double eval_deriv1(int order, double* d, int j)
	{
		double* Hr = Gr(order, j);
		int n = ShapeFunctions(order);
		double s = 0;
		for (int i = 0; i<n; ++i) s += Hr[i] * d[i];
		return s;
	}

	double eval_deriv2(double* d, int j)
	{
		double* Hs = Gs(j);
		int n = Nodes();
		double s = 0;
		for (int i=0; i<n; ++i) s +=  Hs[i]*d[i];
		return s;
	}

	vec3d eval_deriv1(vec3d* d, int j)
	{
		double* Hr = Gr(j);
		int n = Nodes();
		vec3d v(0,0,0);
		for (int i = 0; i<n; ++i) v += d[i]*Hr[i];
		return v;
	}

	vec3d eval_deriv2(vec3d* d, int j)
	{
		double* Hs = Gs(j);
		int n = Nodes();
		vec3d v(0,0,0);
		for (int i = 0; i<n; ++i) v += d[i]*Hs[i];
		return v;
	}

	double eval_deriv1(double* d, double r, double s)
	{
		double Hr[FEElement::MAX_NODES], Hs[FEElement::MAX_NODES];
		shape_deriv(Hr, Hs, r, s);
		int n = Nodes();
		double a = 0;
		for (int i = 0; i<n; ++i) a += Hr[i] * d[i];
		return a;
	}

	double eval_deriv1(int order, double* d, double r, double s)
	{
		double Hr[FEElement::MAX_NODES], Hs[FEElement::MAX_NODES];
		shape_deriv(order, Hr, Hs, r, s);
		int n = ShapeFunctions(order);
		double a = 0;
		for (int i=0; i<n; ++i) a +=  Hr[i]*d[i];
		return a;
	}

	double eval_deriv2(double* d, double r, double s)
	{
		double Hr[FEElement::MAX_NODES], Hs[FEElement::MAX_NODES];
		shape_deriv(Hr, Hs, r, s);
		int n = Nodes();
		double a = 0;
		for (int i=0; i<n; ++i) a +=  Hs[i]*d[i];
		return a;
	}

	double eval_deriv2(int order, double* d, double r, double s)
	{
		double Hr[FEElement::MAX_NODES], Hs[FEElement::MAX_NODES];
		shape_deriv(order, Hr, Hs, r, s);
		int n = ShapeFunctions(order);
		double a = 0;
		for (int i = 0; i<n; ++i) a += Hs[i] * d[i];
		return a;
	}

	void shape_fnc(double* H, double r, double s)
	{
		((FESurfaceElementTraits*)m_pT)->shape_fnc(H, r, s);
	}

	void shape_deriv(double* Gr, double* Gs, double r, double s)
	{
		((FESurfaceElementTraits*)m_pT)->shape_deriv(Gr, Gs, r, s);
	}

	void shape_deriv(int order, double* Gr, double* Gs, double r, double s)
	{
		((FESurfaceElementTraits*)m_pT)->shape_deriv(order, Gr, Gs, r, s);
	}

	void shape_deriv2(double* Grr, double* Grs, double* Gss, double r, double s)
	{
		((FESurfaceElementTraits*)m_pT)->shape_deriv2(Grr, Grs, Gss, r, s);
	}

	void shape_fnc(int order, double* H, double r, double s)
	{
		((FESurfaceElementTraits*)m_pT)->shape_fnc(order, H, r, s);
	}

    //! return number of edges
    int facet_edges() const;
    
    //! return node list of edge
    void facet_edge(int j, int* en) const;

	//! serialize
	void Serialize(DumpStream& ar) override;

	double* Gr(int order, int n) const;
	double* Gs(int order, int n) const;
    
public:
    //! local ID of surface element
	int		m_lid;

	// indices of solid or shell element this surface is a face of
	// For solids, a surface element can be connected to two elements 
	// if the surface is an inside surface. For boundary surfaces
	// the second element index is -1. 
	FEElement*		m_elem[2];
};

