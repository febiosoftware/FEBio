#include "stdafx.h"
#include "Integrate.h"
#include "FESolidDomain.h"

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
