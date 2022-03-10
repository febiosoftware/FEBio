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
#include "FESpringMaterial.h"
#include <FECore/FEElement.h>

vec3d FESpringMaterial::Force(FEDiscreteMaterialPoint& mp)
{
	vec3d e = mp.m_drt; e.unit();

	// calculate spring lengths
	double L0 = mp.m_dr0.norm();
	double Lt = mp.m_drt.norm();
	double DL = Lt - L0;

	// evaluate the spring force
	return e*force(DL);
}

mat3d FESpringMaterial::Stiffness(FEDiscreteMaterialPoint& mp)
{
	vec3d e = mp.m_drt; e.unit();

	// calculate spring lengths
	double L0 = mp.m_dr0.norm();
	double Lt = mp.m_drt.norm();
	double DL = Lt - L0;

	// evaluate the stiffness
	double F = force(DL);
	double E = stiffness(DL);

	if (Lt == 0) { F = 0; Lt = 1; e = vec3d(1, 1, 1); }

	mat3d A; A.zero();
	A[0][0] = ((E - F / Lt)*e.x*e.x + F / Lt);
	A[1][1] = ((E - F / Lt)*e.y*e.y + F / Lt);
	A[2][2] = ((E - F / Lt)*e.z*e.z + F / Lt);

	A[0][1] = A[1][0] = (E - F / Lt)*e.x*e.y;
	A[1][2] = A[2][1] = (E - F / Lt)*e.y*e.z;
	A[0][2] = A[2][0] = (E - F / Lt)*e.x*e.z;

	return A;
}

double FESpringMaterial::StrainEnergy(FEDiscreteMaterialPoint& mp)
{
	// calculate spring lengths
	double L0 = mp.m_dr0.norm();
	double Lt = mp.m_drt.norm();
	double DL = Lt - L0;

	// evaluate the spring force
	return strainEnergy(DL);
}

//-----------------------------------------------------------------------------
// FELinearSpring
//-----------------------------------------------------------------------------

// define the material parameters
BEGIN_FECORE_CLASS(FELinearSpring, FESpringMaterial)
	ADD_PARAMETER(m_E, FE_RANGE_GREATER(0.0), "E");
END_FECORE_CLASS();

double FELinearSpring::force(double dl)
{
	return m_E*dl;
}

double FELinearSpring::stiffness(double dl)
{
	return m_E;
}

double FELinearSpring::strainEnergy(double dl)
{
	return m_E*dl*dl/2;
}

//-----------------------------------------------------------------------------
// FETensionOnlyLinearSpring
//-----------------------------------------------------------------------------

// define the material parameters
BEGIN_FECORE_CLASS(FETensionOnlyLinearSpring, FESpringMaterial)
	ADD_PARAMETER(m_E, FE_RANGE_GREATER(0.0), "E");
END_FECORE_CLASS();

double FETensionOnlyLinearSpring::force(double dl)
{
	if (dl >= 0) return m_E*dl; else return 0;
}

double FETensionOnlyLinearSpring::stiffness(double dl)
{
	return (dl >= 0 ? m_E : 0);
}

double FETensionOnlyLinearSpring::strainEnergy(double dl)
{
	if(dl >= 0) return m_E*dl*dl/2; else return 0;
}

//-----------------------------------------------------------------------------
// FEExperimentalSpring
//-----------------------------------------------------------------------------

// define the material parameters
BEGIN_FECORE_CLASS(FEExperimentalSpring, FESpringMaterial)
	ADD_PARAMETER(m_E , "E");
	ADD_PARAMETER(m_sM, "sM");
	ADD_PARAMETER(m_sm, "sm");
END_FECORE_CLASS();

FEExperimentalSpring::FEExperimentalSpring(FEModel* pfem) : FESpringMaterial(pfem)
{
	m_E = 0.0;
	m_sM = 0.0;
	m_sm = 0.0;
}

double FEExperimentalSpring::force(double dl)
{
	if (dl >= 0.0)
		return m_sM*(1.0 - exp(-m_E*dl / m_sM));
	else
		return -m_sm*(1.0 - exp(m_E*dl / m_sm));
}

double FEExperimentalSpring::stiffness(double dl)
{
	if (dl >= 0.0)
		return m_E*exp(-m_E*dl / m_sM);
	else
		return m_E*exp(m_E*dl / m_sm);
}

double FEExperimentalSpring::strainEnergy(double dl)
{
	if (dl >= 0.0)
		return m_sM*(m_sM/m_E*(exp(-m_E*dl/m_sM) - 1) + dl);
	else
		return m_sm*(m_sm/m_E*(exp(-m_E*dl/m_sm) - 1) - dl);
}
