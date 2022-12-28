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
//!  This class defines a solid element

class FECORE_API FESolidElement : public FEElement
{
public:
	//! default constructor
	FESolidElement() {}

	//! copy constructor
	FESolidElement(const FESolidElement& el);

	//! assignment operator
	FESolidElement& operator = (const FESolidElement& el);

	//! set the element traits
	void SetTraits(FEElementTraits* pt) override;

	double gr(int n) const { return ((FESolidElementTraits*)(m_pT))->gr[n]; }	// integration point coordinate r
	double gs(int n) const { return ((FESolidElementTraits*)(m_pT))->gs[n]; }	// integration point coordinate s
	double gt(int n) const { return ((FESolidElementTraits*)(m_pT))->gt[n]; }	// integration point coordinate t

	double* GaussWeights() const { return &((FESolidElementTraits*)(m_pT))->gw[0]; }			// weights of integration points

	double* Gr(int n) const { return ((FESolidElementTraits*)(m_pT))->m_Gr[n]; }	// shape function derivative to r
	double* Gs(int n) const { return ((FESolidElementTraits*)(m_pT))->m_Gs[n]; }	// shape function derivative to s
	double* Gt(int n) const { return ((FESolidElementTraits*)(m_pT))->m_Gt[n]; }	// shape function derivative to t

	double* Grr(int n) const { return ((FESolidElementTraits*)(m_pT))->Grr[n]; }	// shape function 2nd derivative to rr
	double* Gsr(int n) const { return ((FESolidElementTraits*)(m_pT))->Gsr[n]; }	// shape function 2nd derivative to sr
	double* Gtr(int n) const { return ((FESolidElementTraits*)(m_pT))->Gtr[n]; }	// shape function 2nd derivative to tr

	double* Grs(int n) const { return ((FESolidElementTraits*)(m_pT))->Grs[n]; }	// shape function 2nd derivative to rs
	double* Gss(int n) const { return ((FESolidElementTraits*)(m_pT))->Gss[n]; }	// shape function 2nd derivative to ss
	double* Gts(int n) const { return ((FESolidElementTraits*)(m_pT))->Gts[n]; }	// shape function 2nd derivative to ts

	double* Grt(int n) const { return ((FESolidElementTraits*)(m_pT))->Grt[n]; }	// shape function 2nd derivative to rt
	double* Gst(int n) const { return ((FESolidElementTraits*)(m_pT))->Gst[n]; }	// shape function 2nd derivative to st
	double* Gtt(int n) const { return ((FESolidElementTraits*)(m_pT))->Gtt[n]; }	// shape function 2nd derivative to tt

																					//! values of shape functions
	void shape_fnc(double* H, double r, double s, double t) const { ((FESolidElementTraits*)(m_pT))->shape_fnc(H, r, s, t); }

	//! values of shape function derivatives
	void shape_deriv(double* Hr, double* Hs, double* Ht, double r, double s, double t) const { ((FESolidElementTraits*)(m_pT))->shape_deriv(Hr, Hs, Ht, r, s, t); }

	//! values of shape function second derivatives
	void shape_deriv2(double* Hrr, double* Hss, double* Htt, double* Hrs, double* Hst, double* Hrt, double r, double s, double t) const { ((FESolidElementTraits*)(m_pT))->shape_deriv2(Hrr, Hss, Htt, Hrs, Hst, Hrt, r, s, t); }

	vec3d evaluate(vec3d* v, double r, double s, double t) const;
	double evaluate(double* v, double r, double s, double t) const;

	double* Gr(int order, int n) const;
	double* Gs(int order, int n) const;
	double* Gt(int order, int n) const;

	void Serialize(DumpStream& ar) override;

public:
	std::vector<bool>    m_bitfc;    //!< flag for interface nodes
	std::vector<mat3d>	m_J0i;		//!< inverse of reference Jacobian
};
