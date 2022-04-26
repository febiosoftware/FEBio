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

#include "FEBirzle.h"
#include <math.h>

BEGIN_FECORE_CLASS(FEBirzle, FEElasticMaterial)
	ADD_PARAMETER(m_E, FE_RANGE_GREATER(0.0), "E");
	ADD_PARAMETER(m_v, FE_RANGE_GREATER(0.0), "v");
	ADD_PARAMETER(m_c1, FE_RANGE_GREATER(0.0), "c1");
	ADD_PARAMETER(m_c3, FE_RANGE_GREATER(0.0), "c3");
	ADD_PARAMETER(m_d1, FE_RANGE_GREATER(0.0), "d1");
	ADD_PARAMETER(m_d3, FE_RANGE_GREATER(0.0), "d3");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEBirzle::FEBirzle(FEModel* pfem) : FEElasticMaterial(pfem)
{
}

bool FEBirzle::Init()
{
	if (FEElasticMaterial::Init() == false) return false;

	m_c = m_E/(4*(1+m_v));
	m_b = m_v/(1-2*m_v);

	return true;
}

double FEBirzle::StrainEnergyDensity(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	mat3ds B = pt.LeftCauchyGreen();

	double I1 = B.tr();
	double I3 = B.det();

	return m_c*(I1-3)
		+m_c/m_b*(pow(I3, -m_b)-1)
		+m_c1*pow((pow(I3, -1.0/3.0)*I1-3), m_d1)
		+m_c3*pow((pow(I3, 1.0/3.0)-1), m_d3);
}

mat3ds FEBirzle::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	mat3ds B = pt.LeftCauchyGreen();

	double I1 = B.tr();
	double I3 = B.det();

	double J = pt.m_J;

	mat3dd I(1.0);

	mat3ds s = (2*m_c*B 
		- 2*m_c*pow(I3, -m_b)*I +
		+2*m_c1*m_d1*pow(I3, -1.0/3.0)*pow(pow(I3, -1.0/3.0)*I1-3, m_d1-1)*(B-I1/3.0*I)
		+2.0/3.0*m_c3*m_d3*pow(pow(I3, 1.0/3.0)-1, m_d3-1)*pow(I3,1.0/3.0)*I)/J;

	return s;
}

tens4ds FEBirzle::Tangent(FEMaterialPoint& mp)
{
	// As in the Stress function, we need the data from the FEElasticMaterialPoint
	// class to calculate the tangent.
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// Get the deformation gradient and its determinant
	double J = pt.m_J;

	mat3ds B = pt.LeftCauchyGreen();

	double I1 = B.tr();
	double I3 = B.det();

	// define identity tensor and some useful
	// dyadic products of the identity tensor.
	mat3dd I(1.0);
	tens4ds IxI = dyad1s(I);
	tens4ds I4 = dyad4s(I);
	tens4ds bxI = dyad1s(B,I);
	

	double term = m_c1*m_d1*pow(I3,-1.0/3.0)*pow(pow(I3,-1.0/3.0)*I1-3, m_d1-1);

	tens4ds c = (4*m_c*pow(I3, -m_b)*(m_b*IxI+I4)
		-4.0/3.0*term*(bxI-I1/3.0*IxI)
		+4*m_c1*m_d1*pow(I3,-2.0/3.0)*(m_d1-1)*pow(pow(I3,-1.0/3.0)*I1-3, m_d1-2)*(dyad1s(B)-I1/3.0*bxI+1.0/9.0*I1*I1*IxI)
		+4.0/3.0*term*I1*I4
		+4.0/9.0*m_c3*m_d3*((m_d3-1)*pow(I3, 2.0/3.0)*pow(pow(I3,1.0/3.0)-1, m_d3-2)+pow(pow(I3,1.0/3.0)-1, m_d3-1)*pow(I3,1.0/3.0))*IxI
		-4.0/3.0*m_c3*m_d3*pow(pow(I3,1.0/3.0)-1,m_d3-1)*pow(I3,1.0/3.0)*I4)/J;

	return c;
}
