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
#pragma once
#include "FEUncoupledMaterial.h"

class FEGenericHyperelasticUC : public FEUncoupledMaterial
{
public:
	FEGenericHyperelasticUC(FEModel* fem);

	// initialization
	bool Init() override;

	// serialization
	void Serialize(DumpStream& ar) override;

	//! Deviatoric Cauchy stress
	mat3ds DevStress(FEMaterialPoint& mp) override;

	//! Deviatoric spatial Tangent
	tens4ds DevTangent(FEMaterialPoint& mp) override;

	//! Deviatoric strain energy density
	double DevStrainEnergyDensity(FEMaterialPoint& mp) override;

private:
	bool BuildMathExpressions();

private:
	std::string			m_exp;			// the string with the strain energy expression
	bool				m_printDerivs;	// option to print out derivatives

private:
	MSimpleExpression	m_W;	// deviatoric strain-energy function
	vector<double*>		m_param;	// user parameters

	// strain energy derivatives
	MSimpleExpression	m_W1;
	MSimpleExpression	m_W2;
	MSimpleExpression	m_W11;
	MSimpleExpression	m_W12;
	MSimpleExpression	m_W22;

	DECLARE_FECORE_CLASS();
};
