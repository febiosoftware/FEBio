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
#include "FEUncoupledActiveFiberContraction.h"
#include "FEElasticMaterial.h"
#include <FECore/log.h>

BEGIN_FECORE_CLASS(FEUncoupledActiveFiberContraction, FEActiveContractionMaterial);
	ADD_PARAMETER(m_ascl, "ascl");
	ADD_PARAMETER(m_Tmax, "Tmax");
	ADD_PARAMETER(m_ca0, "ca0");
	ADD_PARAMETER(m_camax, "camax");
	ADD_PARAMETER(m_beta, "beta");
	ADD_PARAMETER(m_l0, "l0");
	ADD_PARAMETER(m_refl, "refl");
END_FECORE_CLASS();

FEUncoupledActiveFiberContraction::FEUncoupledActiveFiberContraction(FEModel* pfem) : FEActiveContractionMaterial(pfem)
{
	m_ascl = 0;
	m_Tmax = 1.0;
	m_ca0 = 1.0;
	m_camax = 0.0;
	m_beta = 1;
	m_l0 = 1;
	m_refl = 1;
}

bool FEUncoupledActiveFiberContraction::Init()
{
	if (FEActiveContractionMaterial::Init() == false) return false;

	// for backward compatibility we set m_camax to m_ca0 if it is not defined
	if (m_camax == 0.0) m_camax = m_ca0;
	if (m_camax <= 0.0) { feLogError("camax must be larger than zero"); return false; }

	return true;
}

mat3ds FEUncoupledActiveFiberContraction::ActiveStress(FEMaterialPoint& mp, const vec3d& a0)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// get the deformation gradient
	mat3d F = pt.m_F;
	double J = pt.m_J;
	double Jm13 = pow(J, -1.0 / 3.0);

	// calculate the current material axis lam*a = F*a0;
	vec3d a = F * a0;

	// normalize material axis and store fiber stretch
	double lam, lamd;
	lam = a.unit();
	lamd = lam * Jm13; // i.e. lambda tilde

	// calculate dyad of a: AxA = (a x a)
	mat3ds AxA = dyad(a);

	// get the activation
	double s = 0.0;
	if (m_ascl > 0)
	{
		// current sarcomere length
		double strl = m_refl * lamd;

		// sarcomere length change
		double dl = strl - m_l0;

		if (dl >= 0)
		{
			// calcium sensitivity
			double eca50i = (exp(m_beta * dl) - 1);

			// ratio of Camax/Ca0
			double rca = m_camax / m_ca0;

			// active fiber stress
			s = m_Tmax * (eca50i / (eca50i + rca * rca)) * m_ascl;
		}
	}

	mat3ds sa = (AxA * (s / J));

	return sa.dev();
}

tens4ds FEUncoupledActiveFiberContraction::ActiveStiffness(FEMaterialPoint& mp, const vec3d& a0)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// get the deformation gradient
	mat3d F = pt.m_F;
	double J = pt.m_J;
	double Jm13 = pow(J, -1.0 / 3.0);

	// calculate the current material axis lam*a = F*a0;
	vec3d a = F * a0;

	// normalize material axis and store fiber stretch
	double lam, lamd;
	lam = a.unit();
	lamd = lam * Jm13; // i.e. lambda tilde
	double lamd2 = lamd * lamd;
	double lamd3 = lamd2 * lamd;

	// calculate dyad of a: AxA = (a x a)
	mat3dd I(1.0);
	mat3ds A = dyad(a);
	tens4ds AxA = dyad1s(A);
	tens4ds IxI = dyad1s(I);

	tens4ds Ax1 = dyad1s(A, I);
	tens4ds i4 = dyad4s(I);

	double sa = 0;
	double dsa = 0;
	if (m_ascl > 0)
	{
		// current sarcomere length
		double strl = m_refl * lamd;

		// sarcomere length change
		double dl = strl - m_l0;

		if (dl >= 0)
		{
			// calcium sensitivity
			double eca50i = (exp(m_beta * dl) - 1);

			// ratio of Camax/Ca0
			double rca = m_camax / m_ca0;
			double r2 = rca * rca;

			// active fiber stress
			double D = eca50i + r2;
			sa = m_Tmax * (eca50i / D) * m_ascl;

			dsa = -2.0 * sa / lamd3 + (1.0 / lamd2) * m_Tmax * m_ascl * m_beta * m_refl * r2 * (eca50i + 1.0) / (D * D);
		}
	}

	tens4ds Jc = Ax1 * (-(2.0 / 3.0) * sa)
		+ (i4 + IxI / 3.0) * (2.0 * sa / 3.0)
		+ (AxA - Ax1 / 3.0 + IxI / 9.0) * (lamd3 * dsa);

	return Jc / J;
}
