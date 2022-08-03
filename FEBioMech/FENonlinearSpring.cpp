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
#include "FENonlinearSpring.h"

BEGIN_FECORE_CLASS(FENonlinearSpringMaterial, FEDiscreteElasticMaterial)
	ADD_PARAMETER(m_scale, "scale");
	ADD_PARAMETER(m_measure, "measure", 0, "elongation\0strain\0stretch\0");
	ADD_PROPERTY(m_F, "force");
END_FECORE_CLASS();


FENonlinearSpringMaterial::FENonlinearSpringMaterial(FEModel* pfem) : FEDiscreteElasticMaterial(pfem)
{
	m_scale = 1.0;
	m_measure = 0;
	m_F = nullptr;
}

vec3d FENonlinearSpringMaterial::Force(FEDiscreteMaterialPoint& mp)
{
	vec3d t = mp.m_drt; t.unit();

	// calculate spring lenghts
	double L0 = mp.m_dr0.norm();
	double Lt = mp.m_drt.norm();
	double DL = Lt - L0;

	// evaluate the deformation measure
	double s = 0.0;
	switch (m_measure)
	{
	case 0: s = DL; break;
	case 1: s = DL / L0; break;
	case 2: s = Lt / L0; break;
	default:
		break;
	}

	// get the force
	double F = m_scale * m_F->value(s);

	// evaluate the spring force
	return t*F;
}

mat3d FENonlinearSpringMaterial::Stiffness(FEDiscreteMaterialPoint& mp)
{
	vec3d e = mp.m_drt; e.unit();

	// calculate spring lengths
	double L0 = mp.m_dr0.norm();
	double Lt = mp.m_drt.norm();
	double DL = Lt - L0;

	// evaluate the deformation measure
	double s = 0.0, ds = 1.0;;
	switch (m_measure)
	{
	case 0: s = DL; break;
	case 1: s = DL / L0; ds = 1.0 / L0; break;
	case 2: s = Lt / L0; ds = 1.0 / L0; break;
	default:
		break;
	}

	// get the force
	double F = m_scale * m_F->value(s);

	// evaluate the stiffness
	double E = m_scale * m_F->derive(s) * ds;

	if (Lt == 0) { F = 0; Lt = 1; e = vec3d(1, 1, 1); }

	mat3d A; A.zero();
	A[0][0] = ((E - F / Lt) * e.x * e.x + F / Lt);
	A[1][1] = ((E - F / Lt) * e.y * e.y + F / Lt);
	A[2][2] = ((E - F / Lt) * e.z * e.z + F / Lt);

	A[0][1] = A[1][0] = (E - F / Lt) * e.x * e.y;
	A[1][2] = A[2][1] = (E - F / Lt) * e.y * e.z;
	A[0][2] = A[2][0] = (E - F / Lt) * e.x * e.z;

	return A;
}

double FENonlinearSpringMaterial::StrainEnergy(FEDiscreteMaterialPoint& mp)
{
	// TODO: implement this
	assert(false);
	return 0.0;
}
