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
#include "FEUncoupledActiveContraction.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEUncoupledActiveContraction, FEUncoupledMaterial)
	ADD_PARAMETER(m_Tmax , "Tmax" );
	ADD_PARAMETER(m_ca0  , "ca0"  );
	ADD_PARAMETER(m_camax, "camax");
	ADD_PARAMETER(m_beta , "beta" );
	ADD_PARAMETER(m_l0   , "l0"   );
	ADD_PARAMETER(m_refl , "refl" );

	ADD_PROPERTY(m_Q, "mat_axis")->SetFlags(FEProperty::Optional);

END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEUncoupledActiveContraction::FEUncoupledActiveContraction(FEModel* pfem) : FEUncoupledMaterial(pfem)
{

}

//-----------------------------------------------------------------------------
mat3ds FEUncoupledActiveContraction::DevStress(FEMaterialPoint &mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// deformation gradient
	mat3d& F = pt.m_F;
	double J = pt.m_J;
	double Jm13 = pow(J, -1.0/3.0);

	// get the local coordinate systems
	mat3d Q = GetLocalCS(mp);

	// get the initial fiber direction
	vec3d a0 = Q.col(0);

	// calculate the current material axis lam*a = F*a0;
	vec3d a = F*a0;

	// normalize material axis and store fiber stretch
	double lam = a.unit();
	double lamd = lam*Jm13; // i.e. lambda tilde

	// current sarcomere length
	double strl = m_refl*lamd;

	// sarcomere length change
	double dl = strl - m_l0;

	// calculate stress
	mat3ds s; s.zero();
	if (dl >= 0)
	{
		// calcium sensitivity
		double eca50i = (exp(m_beta*dl) - 1);

		// ratio of Camax/Ca0
		double rca = m_camax/m_ca0;

		// active fiber stress
		double saf = m_Tmax*(eca50i / ( eca50i + rca*rca ));

		// calculate dyad of a: AxA = (a x a)
		mat3ds AxA = dyad(a);

		// add saf*(a x a)
		s += AxA*saf;
	}

	return s;
}

//-----------------------------------------------------------------------------
// \todo Implement this
tens4ds FEUncoupledActiveContraction::DevTangent(FEMaterialPoint &pt)
{
	return tens4ds(0.0);
}
