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
//!  This class defines the shell element. 

//! A shell element is similar to a surface
//! element except that it has a thickness. 

class FECORE_API FEShellElement : public FEElement
{
public:
	FEShellElement();

	//! copy constructor
	FEShellElement(const FEShellElement& el);

	//! assignment operator
	FEShellElement& operator = (const FEShellElement& el);

	virtual void SetTraits(FEElementTraits* ptraits) override;

	double gr(int n) { return ((FEShellElementTraits*)(m_pT))->gr[n]; }
	double gs(int n) { return ((FEShellElementTraits*)(m_pT))->gs[n]; }
	double gt(int n) { return ((FEShellElementTraits*)(m_pT))->gt[n]; }

	double* GaussWeights() { return &((FEShellElementTraits*)(m_pT))->gw[0]; }	// weights of integration points

	double* Hr(int n) { return ((FEShellElementTraits*)(m_pT))->Hr[n]; }	// shape function derivative to r
	double* Hs(int n) { return ((FEShellElementTraits*)(m_pT))->Hs[n]; }	// shape function derivative to s

																			//! values of shape functions
	void shape_fnc(double* H, double r, double s) const { ((FEShellElementTraits*)(m_pT))->shape_fnc(H, r, s); }

	//! values of shape function derivatives
	void shape_deriv(double* Hr, double* Hs, double r, double s) const { ((FEShellElementTraits*)(m_pT))->shape_deriv(Hr, Hs, r, s); }

	//! serialize data associated with this element
	void Serialize(DumpStream &ar) override;

public:
	std::vector<double>	m_h0;	//!< initial shell thicknesses
	std::vector<double>	m_ht;	//!< current shell thickness
	std::vector<vec3d>	m_d0;   //!< initial shell director

	std::vector<vec3d>	m_g0[3];//!< reference covariant base vectors
	std::vector<vec3d>	m_gt[3];//!< current covariant base vectors
	std::vector<vec3d>	m_gp[3];//!< previous covariant base vectors

	std::vector<vec3d>	m_G0[3];//!< reference contravariant base vectors
	std::vector<vec3d>	m_Gt[3];//!< current contravariant base vectors

	// indices of solid elements this shell element is attached to.
	// the first element is attached to the back of the shell
	// and the second element is attached to the front.
	// the index is -1 if no solid is attached on that side.
	int        m_elem[2];
};

//-----------------------------------------------------------------------------
// Shell element used by old shell formulation
class FECORE_API FEShellElementOld : public FEShellElement
{
public:
	FEShellElementOld();

	//! copy constructor
	FEShellElementOld(const FEShellElementOld& el);

	//! assignment operator
	FEShellElementOld& operator = (const FEShellElementOld& el);

	// set the element traits class
	void SetTraits(FEElementTraits* ptraits) override;

	//! serialize data associated with this element
	void Serialize(DumpStream &ar) override;

public:
	std::vector<vec3d>	m_D0;	//!< initial shell directors
};

//-----------------------------------------------------------------------------
// Shell element used by new shell formulations
class FECORE_API FEShellElementNew : public FEShellElement
{
public:
	FEShellElementNew();

	//! copy constructor
	FEShellElementNew(const FEShellElementNew& el);

	//! assignment operator
	FEShellElementNew& operator = (const FEShellElementNew& el);

	// set the element traits class
	void SetTraits(FEElementTraits* ptraits) override;

	//! serialize data associated with this element
	void Serialize(DumpStream &ar) override;

public: // EAS parameters

	matrix          m_Kaai;
	matrix          m_fa;
	matrix          m_alpha;
	matrix          m_alphat;
	matrix          m_alphai;
	std::vector<matrix>  m_Kua;
	std::vector<matrix>  m_Kwa;
	std::vector<mat3ds>  m_E;
};

