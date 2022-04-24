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
#include "FEActiveFiberStressUC.h"
#include "FEElasticMaterial.h"
#include <FECore/log.h>

//=====================================================================================

BEGIN_FECORE_CLASS(FEActiveFiberStressUC, FEUncoupledMaterial);
	ADD_PARAMETER(m_smax, "smax")->setUnits(UNIT_PRESSURE);
	ADD_PARAMETER(m_ac, "activation");

	ADD_PROPERTY(m_stl, "stl", FEProperty::Optional);
	ADD_PROPERTY(m_stv, "stv", FEProperty::Optional);

	ADD_PROPERTY(m_Q, "mat_axis")->SetFlags(FEProperty::Optional);

END_FECORE_CLASS();

FEActiveFiberStressUC::FEActiveFiberStressUC(FEModel* fem) : FEUncoupledMaterial(fem)
{
	m_smax = 0.0;
	m_ac = 0.0;

	m_stl = nullptr;
	m_stv = nullptr;
}

mat3ds FEActiveFiberStressUC::DevStress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	double dt = 0.0;

	mat3d& F = pt.m_F;
	double J = pt.m_J;
	double J13 = pow(J, 1.0 / 3.0);

	mat3d Q = GetLocalCS(mp);

	vec3d a0 = Q.col(0);

	vec3d a = F*a0;
	double lam = a.unit() / J13;

	double stl = (m_stl ? m_stl->value(lam) : 1.0);
	double v = 0;// (lam - lamp) / dt;
	double stv = (m_stv ? m_stl->value(v) : 1.0);

	mat3ds A = dyad(a);

	double sf = m_ac*m_smax*stl*stv;

	return A.dev()*(sf / J);
}

tens4ds FEActiveFiberStressUC::DevTangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	double dt = 1.0;

	mat3d& F = pt.m_F;
	double J = pt.m_J;
	double J13 = pow(J, 1.0 / 3.0);

	mat3d Q = GetLocalCS(mp);

	vec3d a0 = Q.col(0);

	vec3d a = F*a0;
	double lam = a.unit() / J13;

	mat3ds A = dyad(a);

	mat3dd I(1.0);
	tens4ds IxI = dyad1s(I);
	tens4ds I4 = dyad4s(I);

	tens4ds AA = dyad1s(A.dev());
	tens4ds AI = dyad1s(A, I);

	double stl = (m_stl ? m_stl->value(lam) : 1.0);
	double v = 0;// (lam - lamp) / dt;
	double stv = (m_stv ? m_stl->value(v) : 1.0);

	double dstl = (m_stl ? m_stl->derive(lam) : 0.0);
	double dstv = (m_stv ? m_stv->derive(v) / dt : 0.0);

	double sf = m_ac*m_smax*stl*stv;
	double sf_l = m_ac*m_smax*(dstl*stv + stl*dstv);

	tens4ds cff = AA*(lam*sf_l - 2.0*sf);
	tens4ds cf = (AI - I4 - IxI/3.0)*(2.0*sf/3.0);

	tens4ds c = (cff - cf) / J;

	return c;
}
