#include "stdafx.h"
#include "Integrate.h"
#include "FESolidDomain.h"

void IntegrateBDB(FESolidDomain& dom, FESolidElement& el, double D, matrix& ke)
{
	int ne = el.Nodes();
	int ni = el.GaussPoints();

	const int EN = FEElement::MAX_NODES;
	double Gx[EN], Gy[EN], Gz[EN];
	double Ji[3][3];
	double Gi[3], Gj[3];
	double DB[3];

	// zero stiffness matrix
	ke.zero();

	// loop over all integration points
	const double *gw = el.GaussWeights();
	for (int n = 0; n<ni; ++n)
	{
		// calculate jacobian
		double detJt = dom.invjact(el, Ji, n);

		// evaluate the conductivity
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);

		for (int i = 0; i<ne; ++i)
		{
			double Gr = el.Gr(n)[i];
			double Gs = el.Gs(n)[i];
			double Gt = el.Gt(n)[i];

			// calculate global gradient of shape functions
			// note that we need the transposed of Ji, not Ji itself !
			Gx[i] = Ji[0][0] * Gr + Ji[1][0] * Gs + Ji[2][0] * Gt;
			Gy[i] = Ji[0][1] * Gr + Ji[1][1] * Gs + Ji[2][1] * Gt;
			Gz[i] = Ji[0][2] * Gr + Ji[1][2] * Gs + Ji[2][2] * Gt;
		}

		for (int i = 0; i<ne; ++i)
		{
			Gi[0] = Gx[i];
			Gi[1] = Gy[i];
			Gi[2] = Gz[i];

			for (int j = 0; j<ne; ++j)
			{
				Gj[0] = Gx[j];
				Gj[1] = Gy[j];
				Gj[2] = Gz[j];

				DB[0] = D*Gj[0];
				DB[1] = D*Gj[1];
				DB[2] = D*Gj[2];

				ke[i][j] += (Gi[0] * DB[0] + Gi[1] * DB[1] + Gi[2] * DB[2])*detJt*gw[n];
			}
		}
	}
}
