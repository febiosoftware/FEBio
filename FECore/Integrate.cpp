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
#include "Integrate.h"
#include "FESolidDomain.h"
#include "FELinearSystem.h"
#include "FESolver.h"

//-----------------------------------------------------------------------------
void FECORE_API IntegrateBDB(FESolidDomain& dom, FESolidElement& el, double D, matrix& ke)
{
	// vector to store global shape functions
	const int EN = FEElement::MAX_NODES;
	vec3d G[EN];

	// loop over all integration points
	const double *gw = el.GaussWeights();
	int ne = el.Nodes();
	int ni = el.GaussPoints();
	for (int n = 0; n<ni; ++n)
	{
		// calculate jacobian
		double detJt = dom.ShapeGradient(el, n, G);

		// form the matrix
		for (int i = 0; i<ne; ++i)
		{
			for (int j = 0; j<ne; ++j)
			{
				ke[i][j] += (G[i]*G[j])*(D*detJt*gw[n]);
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FECORE_API IntegrateBDB(FESolidDomain& dom, FESolidElement& el, const mat3ds& D, matrix& ke)
{
	// vector to store global shape functions
	const int EN = FEElement::MAX_NODES;
	vec3d G[EN];

	// loop over all integration points
	const double *gw = el.GaussWeights();
	int ne = el.Nodes();
	int ni = el.GaussPoints();
	for (int n = 0; n<ni; ++n)
	{
		// calculate jacobian
		double detJt = dom.ShapeGradient(el, n, G);

		// form the matrix
		for (int i = 0; i<ne; ++i)
		{
			for (int j = 0; j<ne; ++j)
			{
				ke[i][j] += (G[i] * (D * G[j]))*(detJt*gw[n]);
			}
		}
	}
}

//-----------------------------------------------------------------------------
FECORE_API void IntegrateBDB(FESolidDomain& dom, FESolidElement& el, std::function<mat3ds(const FEMaterialPoint& mp)> D, matrix& ke)
{
	// vector to store global shape functions
	const int EN = FEElement::MAX_NODES;
	vec3d G[EN];

	// loop over all integration points
	const double *gw = el.GaussWeights();
	int ne = el.Nodes();
	int ni = el.GaussPoints();
	for (int n = 0; n<ni; ++n)
	{
		// get the material point
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);

		// calculate jacobian
		double detJt = dom.ShapeGradient(el, n, G);

		// calculate D at this point
		mat3ds Dn = D(mp);

		// form the matrix
		for (int i = 0; i<ne; ++i)
		{
			for (int j = 0; j<ne; ++j)
			{
				ke[i][j] += (G[i] * (Dn * G[j]))*(detJt*gw[n]);
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FECORE_API IntegrateNCN(FESolidDomain& dom, FESolidElement& el, double C, matrix& ke)
{
	// number of nodes
	int ne = el.Nodes();

	// jacobian
	double Ji[3][3];

	// loop over all integration points
	const double *gw = el.GaussWeights();
	int ni = el.GaussPoints();
	for (int n = 0; n<ni; ++n)
	{
		// calculate jacobian
		double detJt = dom.invjact(el, Ji, n);

		// shape function values at integration point n
		double* H = el.H(n);

		for (int i = 0; i<ne; ++i)
		{
			for (int j = 0; j<ne; ++j)
			{
				ke[i][j] += H[i] * H[j]*C*detJt*gw[n];
			}
		}
	}
}

//-----------------------------------------------------------------------------
// Generice integrator class for solid domains
FECORE_API void AssembleSolidDomain(FESolidDomain& dom, FELinearSystem& ls, std::function<void(FESolidElement& el, matrix& ke)> elementIntegrand)
{
	// loop over all elements in domain
	int NE = dom.Elements();
#pragma omp parallel for shared (NE)
	for (int i = 0; i<NE; ++i)
	{
		FESolidElement& el = dom.Element(i);
		int ndofs = dom.GetElementDofs(el);

		// build the element stiffness matrix
		FEElementMatrix ke(ndofs, ndofs);
		elementIntegrand(el, ke);

		// set up the LM matrix
		vector<int> lm;
		dom.UnpackLM(el, lm);

		// assemble into global matrix
		ke.SetNodes(el.m_node);
		ke.SetIndices(lm);
		ls.Assemble(ke);
	}
}

//-----------------------------------------------------------------------------
// Generic integrator class for solid domains
FECORE_API void AssembleSolidDomain(FESolidDomain& dom, FEGlobalVector& R, std::function<void(FESolidElement& el, vector<double>& fe)> elementIntegrand)
{
	int NE = dom.Elements();
	for (int i = 0; i < NE; ++i)
	{
		FESolidElement& el = dom.Element(i);
		int ndofs = dom.GetElementDofs(el);

		// get element contribution
		vector<double> fe(ndofs, 0.0);
		elementIntegrand(el, fe);

		// assemble into RHS
		vector<int> lm;
		dom.UnpackLM(el, lm);
		R.Assemble(lm, fe);
	}
}

//-----------------------------------------------------------------------------
// Generice integrator class for solid domains
FECORE_API void IntegrateSolidDomain(FESolidDomain& dom, FELinearSystem& ls, std::function<void(FEMaterialPoint& mp, matrix& ke)> elementIntegrand)
{
	// loop over all elements in domain
	int NE = dom.Elements();
#pragma omp parallel for shared (NE)
	for (int i = 0; i<NE; ++i)
	{
		FESolidElement& el = dom.Element(i);
		int ndofs = dom.GetElementDofs(el);

		// build the element stiffness matrix
		FEElementMatrix ke(ndofs, ndofs);
		ke.zero();
		matrix kn(ndofs, ndofs);

		// loop over all integration points
		int nint = el.GaussPoints();
		double* w = el.GaussWeights();
		for (int n = 0; n < nint; ++n)
		{
			FEMaterialPoint& mp = *el.GetMaterialPoint(n);

			// evaluate the integration point's contribution
			elementIntegrand(mp, kn);

			// add it to the element matrix
			ke.adds(kn, w[n]);
		}

		// set up the LM matrix
		vector<int> lm;
		dom.UnpackLM(el, lm);

		// assemble into global matrix
		ke.SetNodes(el.m_node);
		ke.SetIndices(lm);
		ls.Assemble(ke);
	}
}