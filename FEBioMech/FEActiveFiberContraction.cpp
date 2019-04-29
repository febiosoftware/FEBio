/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in 
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
#include "FEActiveFiberContraction.h"
#include "FEElasticMaterial.h"
#include <FECore/log.h>

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEActiveFiberContraction, FEMaterial);
	ADD_PARAMETER(m_ascl , "ascl");
	ADD_PARAMETER(m_Tmax , "Tmax");
	ADD_PARAMETER(m_ca0  , "ca0");
	ADD_PARAMETER(m_camax, "camax");
	ADD_PARAMETER(m_beta , "beta");
	ADD_PARAMETER(m_l0   , "l0");
	ADD_PARAMETER(m_refl , "refl");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEActiveFiberContraction::FEActiveFiberContraction(FEModel* pfem) : FEMaterial(pfem)
{
	m_ascl = 0;
	m_Tmax = 1.0;
	m_ca0 = 1.0;
	m_camax = 0.0;
}

//-----------------------------------------------------------------------------
bool FEActiveFiberContraction::Init()
{
	if (FEMaterial::Init() == false) return false;

	// for backward compatibility we set m_camax to m_ca0 if it is not defined
	if (m_camax == 0.0) m_camax = m_ca0;
	if (m_camax <= 0.0) { feLogError("camax must be larger than zero"); return false; }

	return true;
}

//-----------------------------------------------------------------------------
mat3ds FEActiveFiberContraction::FiberStress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// get the deformation gradient
	mat3d F = pt.m_F;
	double J = pt.m_J;
	double Jm13 = pow(J, -1.0 / 3.0);

	// get the local coordinate systems
	mat3d Q = GetLocalCS(mp);

	// get the initial fiber direction
	vec3d a0 = Q.col(0);

	// calculate the current material axis lam*a = F*a0;
	vec3d a = F*a0;

	// normalize material axis and store fiber stretch
	double lam, lamd;
	lam = a.unit();
	lamd = lam*Jm13; // i.e. lambda tilde

	// calculate dyad of a: AxA = (a x a)
	mat3ds AxA = dyad(a);

	// get the activation
	double saf = 0.0;
	if (m_ascl > 0)
	{
		double ctenslm = m_ascl;

		// current sarcomere length
		double strl = m_refl*lamd;

		// sarcomere length change
		double dl = strl - m_l0;

		if (dl >= 0)
		{
			// calcium sensitivity
			double eca50i = (exp(m_beta*dl) - 1);

			// ratio of Camax/Ca0
			double rca = m_camax / m_ca0;

			// active fiber stress
			saf = m_Tmax*(eca50i / (eca50i + rca*rca))*ctenslm;
		}
	}
	return AxA*saf;
}

//-----------------------------------------------------------------------------
tens4ds FEActiveFiberContraction::FiberStiffness(FEMaterialPoint& mp)
{
	/*	if (lcna >= 0)
	{
	double ctenslm = m_plc->Value();

	// current sarcomere length
	double strl = refl*lamd;

	// sarcomere length change
	double dl = strl - l0;

	if (dl >= 0) W44 += J*2*beta*refl*exp(-beta*dl);
	}
	*/
	return tens4ds(0.0);
}
