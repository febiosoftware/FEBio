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
#include "FEDiscreteElementMaterial.h"
#include <FECore/FEModel.h>

BEGIN_FECORE_CLASS(FEDiscreteContractileMaterial, FEDiscreteElasticMaterial)
	ADD_PARAMETER(m_Vmax, "Vmax"  );
	ADD_PARAMETER(m_ac  , "ac"  );
	ADD_PARAMETER(m_Fmax, "Fmax");
	ADD_PARAMETER(m_Ksh , "Ksh" );
	ADD_PARAMETER(m_Lmax, "Lmax");
	ADD_PARAMETER(m_L0  , "L0" );

	ADD_PROPERTY(m_Sv , "Sv" , FEProperty::Optional);
	ADD_PROPERTY(m_Ftl, "Ftl", FEProperty::Optional);
	ADD_PROPERTY(m_Ftv, "Fvl", FEProperty::Optional);
END_FECORE_CLASS();

FEDiscreteContractileMaterial::FEDiscreteContractileMaterial(FEModel* fem) : FEDiscreteElasticMaterial(fem)
{
	m_Vmax = 1.0;
	m_ac   = 0.0;
	m_Fmax = 1.0;
	m_L0   = 0.0;

	m_Sv = nullptr;
	m_Ftl = nullptr;
	m_Ftv = nullptr;
}

double FEDiscreteContractileMaterial::passive_force(double L, double V)
{
	return (L > 1.0 ? m_Fmax * ((exp(m_Ksh*(L - 1.0) / m_Lmax) - 1.0) / (exp(m_Ksh) - 1.0)) : 0.0);
}

double FEDiscreteContractileMaterial::passive_force_deriv_L(double L, double V)
{
	if (L < 1.0) return 0.0;
	double dF = (m_Fmax / (exp(m_Ksh) - 1.0)) * ((m_Ksh / m_Lmax)*(exp(m_Ksh*(L - 1.0) / m_Lmax)));
	return dF;
}

double FEDiscreteContractileMaterial::passive_force_deriv_V(double L, double V)
{
	return 0.0;
}

double FEDiscreteContractileMaterial::active_force(double L, double V)
{
	double Ftl = (m_Ftl ? m_Ftl->value(L) : 1.0);
	double Ftv = (m_Ftv ? m_Ftv->value(V) : 1.0);

	return m_ac*m_Fmax*Ftl*Ftv;
}

double FEDiscreteContractileMaterial::active_force_deriv_L(double L, double V)
{
	double dFtl = (m_Ftl ? m_Ftl->derive(L) : 0.0);
	double Ftv = (m_Ftv ? m_Ftv->value(V) : 1.0);
	return m_ac*m_Fmax*dFtl*Ftv;
}

double FEDiscreteContractileMaterial::active_force_deriv_V(double L, double V)
{
	double Ftl = (m_Ftl ? m_Ftl->value(L) : 1.0);
	double dFtv = (m_Ftv ? m_Ftv->derive(V) : 0.0);
	return m_ac*m_Fmax*Ftl*dFtv;
}

double FEDiscreteContractileMaterial::force_deriv_L(FEDiscreteMaterialPoint& mp)
{
	vec3d e0 = mp.m_dr0; double L0 = e0.unit();
	vec3d et = mp.m_drt; double Lm = et.unit();
	if (m_L0 != 0.0) L0 = m_L0;

	// stretch
	double L = Lm / L0;

	// normalized velocity
	double Vm = mp.m_dvt*et;
	double Sv = (m_Sv ? m_Sv->value(m_ac) : 1.0);
	double V = Vm / (m_Vmax * Sv);

	return passive_force_deriv_L(L, V)/L0 + active_force_deriv_L(L, V)/L0;
}

double FEDiscreteContractileMaterial::force_deriv_V(FEDiscreteMaterialPoint& mp)
{
	vec3d e0 = mp.m_dr0; double L0 = e0.unit();
	vec3d et = mp.m_drt; double Lm = et.unit();
	if (m_L0 != 0.0) L0 = m_L0;

	// stretch
	double L = Lm / L0;

	// normalized velocity
	double Vm = mp.m_dvt*et;
	double Sv = (m_Sv ? m_Sv->value(m_ac) : 1.0);
	double V = Vm / (m_Vmax * Sv);

	return passive_force_deriv_V(L, V) / (m_Vmax * Sv) + active_force_deriv_V(L, V) / (m_Vmax * Sv);
}

double FEDiscreteContractileMaterial::force(FEDiscreteMaterialPoint& mp)
{
	vec3d e0 = mp.m_dr0; double L0 = e0.unit();
	vec3d et = mp.m_drt; double Lm = et.unit();
	if (m_L0 != 0.0) L0 = m_L0;

	// stretch
	double L = Lm / L0;

	// normalized velocity
	double Vm = mp.m_dvt*et;
	double Sv = (m_Sv ? m_Sv->value(m_ac) : 1.0);
	double V = Vm / (m_Vmax * Sv);

	// passive element
	double Fp = passive_force(L, V);

	// active element
	double Fc = active_force(L, V);

	// add them up together
	double F = Fp + Fc;

	// return
	return F;
}

// evaluate the force at a discrete element
vec3d FEDiscreteContractileMaterial::Force(FEDiscreteMaterialPoint& mp)
{
	vec3d e = mp.m_drt; e.unit();
	double F = force(mp);
	return e*F;
}

// evaluate the stiffness at a discrete element (= dF / dr)
mat3d FEDiscreteContractileMaterial::Stiffness(FEDiscreteMaterialPoint& mp)
{
	double F = force(mp);

	vec3d e = mp.m_drt;
	double L = e.unit();

	mat3ds exe = dyad(e);
	mat3dd I(1.0);

	mat3ds E = (I - exe)/L;

	vec3d v = mp.m_dvt;
	double V = mp.m_dvt * e;

	double Fl = force_deriv_L(mp);
	double Fv = force_deriv_V(mp);

	double dt = GetFEModel()->GetTime().timeIncrement;

	vec3d f = e*(-Fl) - (e/dt + v/L + e*(L*V))*Fv;

	return E*F - (e & f);
}
