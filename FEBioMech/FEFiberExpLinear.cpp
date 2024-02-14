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
#include <limits>
#include "FEFiberExpLinear.h"
#include <FECore/expint_Ei.h>

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEFiberExpLinear, FEFiberMaterial)
	ADD_PARAMETER(m_c3  , FE_RANGE_GREATER_OR_EQUAL(0.0), "c3")->setUnits(UNIT_PRESSURE);
	ADD_PARAMETER(m_c4  , FE_RANGE_GREATER_OR_EQUAL(0.0), "c4");
    ADD_PARAMETER(m_c5  , FE_RANGE_GREATER_OR_EQUAL(0.0), "c5")->setUnits(UNIT_PRESSURE);
	ADD_PARAMETER(m_lam1, FE_RANGE_GREATER_OR_EQUAL(1.0), "lambda");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FEFiberExpLinear::FEFiberExpLinear(FEModel* pfem) : FEFiberMaterial(pfem)
{
	m_c3 = 0;
	m_c4 = 0;
    m_c5 = 0;
	m_lam1 = 1.0;
	m_epsf = 0.0;
}

//-----------------------------------------------------------------------------
//! Calculate the fiber stress
mat3ds FEFiberExpLinear::FiberStress(FEMaterialPoint& mp, const vec3d& a0)
{
	// get the material point data
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

    // get the spatial fiber vector and stretch
    vec3d a = pt.m_F*a0;
    double l = a.unit();
    
    // fiber stress
    mat3ds s; s.zero();
    
    if (l >= 1.0)
    {
        // other stuff we need
        mat3ds A = dyad(a);
        double J = pt.m_J;
        
        double c3 = m_c3(mp);
        double c4 = m_c4(mp);
        double c5 = m_c5(mp);
        double lam1 = m_lam1(mp);
        
        if (c3 == 0) {
            c3 = c5/c4*exp(-c4*(lam1-1));
            double c6 = c3*(exp(c4*(lam1-1))*(1-c4*lam1)-1);
            
            // calculate fiber stress
            double sn = 0.0;
            if (l < lam1)
                sn = c3*(exp(c4*(l - 1.0)) - 1.0);
            else
                sn = c5*l + c6;
            s += A*(sn / J);
        }
        else {
            // calculate fiber stress
            double Wl = 0.0;
            if (l < lam1)
            {
                Wl = c3*(exp(c4*(l - 1.0)) - 1.0);
            }
            else
            {
                double c6 = c3*(exp(c4*(lam1 - 1.0)) - 1.0) - c5*lam1;
                Wl = c5*l + c6;
            }
            s += A*(Wl / J);
        }
    }
	return s;
}

//-----------------------------------------------------------------------------
//! Calculate the fiber tangent
tens4ds FEFiberExpLinear::FiberTangent(FEMaterialPoint& mp, const vec3d& a0)
{
	double eps = m_epsf * std::numeric_limits<double>::epsilon();

	// get material point data
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

    // get the spatial fiber axis
    vec3d a = pt.m_F*a0;
    double l = a.unit();

    // fiber tangent
    tens4ds c(0.0);
    
    if (l >= 1.0 + eps)
    {
        // Invariants of B (= invariants of C)
        double J = pt.m_J;
        double I4 = l*l;
        
        // some useful tensors
        mat3dd I(1.0);
        mat3ds A = dyad(a);
        tens4ds AxA = dyad1s(A);
        double c3 = m_c3(mp);
        double c4 = m_c4(mp);
        double c5 = m_c5(mp);
        double lam1 = m_lam1(mp);

        if (c3 == 0) {
            c3 = c5/c4*exp(-c4*(lam1-1));
            double c6 = c3*(exp(c4*(lam1-1))*(1-c4*lam1)-1);
            
            double cn = 0;
            if (l < lam1)
                cn = c3*(2+exp(c4*(l-1))*(c4*l-2));
            else
                cn = -c5*l - 2*c6;
            
            c += AxA*(cn / J);
        }
        else {
            double Fl = 0.0, Fll = 0.0;
            if (l < lam1)
            {
                Fl = c3*(exp(c4*(l - 1.0)) - 1.0) / l;
                Fll = -c3*(exp(c4*(l - 1.0)) - 1.0) / (l*l) + c3*c4*exp(c4*(l - 1.0)) / l;
            }
            else
            {
                double c6 = c3*(exp(c4*(lam1 - 1.0)) - 1.0) - c5*lam1;
                Fl = c5 + c6 / l;
                Fll = -c6 / (l*l);
            }

            double W44 = (Fll - Fl / l) / (4 * l*l);

            c += AxA*(4.0*W44*I4*I4 / J);
        }
    }
	return c;
}

//-----------------------------------------------------------------------------
//! Calculate the fiber strain energy density
double FEFiberExpLinear::FiberStrainEnergyDensity(FEMaterialPoint& mp, const vec3d& a0)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// get the deformation gradient
	mat3d F = pt.m_F;

	// calculate the current material axis lam*a = F*a0;
	vec3d a = F * a0;

	// normalize material axis and store fiber stretch
	double lam = a.unit();
	
	// strain energy density
	double sed = 0.0;
    if (lam >= 1)
    {
        double c3 = m_c3(mp);
        double c4 = m_c4(mp);
        double c5 = m_c5(mp);
        double lam1 = m_lam1(mp);
        
        if (c3 == 0) {
            c3 = c5/c4*exp(-c4*(lam1-1));
            double c6 = c3*(exp(c4*(lam1-1))*(1-c4*lam1)-1);
            if (lam < lam1)
                sed = c3*exp(-c4)*(expint_Ei(c4*lam) - expint_Ei(c4)) - c3*log(lam);
            else
                sed = c5 * (lam - lam1) + c6 * log(lam/lam1)
                + c3*exp(-c4)*(expint_Ei(c4*lam1) - expint_Ei(c4)) - c3*log(lam1);
        }
        else {
            if (lam < lam1)
            {
                sed = c3 * (exp(-c4)*
                    (expint_Ei(c4*lam) - expint_Ei(c4))
                    - log(lam));
            }
            else
            {
                double c6 = c3 * (exp(c4*(lam1 - 1)) - 1) - c5 * lam1;
                sed = c5 * (lam - 1) + c6 * log(lam);
            }
        }
	}

	return sed;
}


//============================================================================================
// define the material parameters
BEGIN_FECORE_CLASS(FEElasticFiberExpLinear, FEElasticFiberMaterial)
    ADD_PARAMETER(m_fib.m_c3  , FE_RANGE_GREATER_OR_EQUAL(0.0), "c3")->setUnits(UNIT_PRESSURE);
    ADD_PARAMETER(m_fib.m_c4  , FE_RANGE_GREATER_OR_EQUAL(0.0), "c4");
    ADD_PARAMETER(m_fib.m_c5  , FE_RANGE_GREATER_OR_EQUAL(0.0), "c5")->setUnits(UNIT_PRESSURE);
    ADD_PARAMETER(m_fib.m_lam1, FE_RANGE_GREATER_OR_EQUAL(1.0), "lambda");
END_FECORE_CLASS();
