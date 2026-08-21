/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2026 University of Utah, The Trustees of Columbia University in
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
#include "FEScriptedElasticMaterial.h"

FEScriptedElasticMaterial::FEScriptedElasticMaterial(FEModel* pfem) : FEScripted<FEElasticMaterial>(pfem)
{
	ScriptContext sc;
	sc.returnType = FEValueType::Mat3d;
	sc.addVariable("C", FEValueType::Mat3d, true, true);
	SetScriptContext(sc);
}

FEScriptedElasticMaterial::~FEScriptedElasticMaterial()
{

}

mat3ds FEScriptedElasticMaterial::Stress(FEMaterialPoint& pt)
{
	FEElasticMaterialPoint& ep = *pt.ExtractData<FEElasticMaterialPoint>();

	std::vector<FEValue> vars(1);
	vars[0] = ep.RightCauchyGreen();

	mat3ds S = Value(pt, vars).toMat3d().sym();
	mat3ds s = ep.push_forward(S);
	return s;
}

tens4ds FEScriptedElasticMaterial::Tangent(FEMaterialPoint& pt)
{
	FEElasticMaterialPoint& ep = *pt.ExtractData<FEElasticMaterialPoint>();

	std::vector<FEValue> vars(2); // we need one more for the directional derivative
	vars[0] = ep.RightCauchyGreen();

	// calculate derivatives
	vars[1] = mat3d(1, 0, 0, 0, 0, 0, 0, 0, 0);	mat3ds dSdC00 = DerivValue(pt, vars, 0).toMat3d().sym();
	vars[1] = mat3d(0, 1, 0, 0, 0, 0, 0, 0, 0);	mat3ds dSdC01 = DerivValue(pt, vars, 0).toMat3d().sym();
	vars[1] = mat3d(0, 0, 1, 0, 0, 0, 0, 0, 0);	mat3ds dSdC02 = DerivValue(pt, vars, 0).toMat3d().sym();
	vars[1] = mat3d(0, 0, 0, 0, 1, 0, 0, 0, 0);	mat3ds dSdC11 = DerivValue(pt, vars, 0).toMat3d().sym();
	vars[1] = mat3d(0, 0, 0, 0, 0, 1, 0, 0, 0);	mat3ds dSdC12 = DerivValue(pt, vars, 0).toMat3d().sym();
	vars[1] = mat3d(0, 0, 0, 0, 0, 0, 0, 0, 1);	mat3ds dSdC22 = DerivValue(pt, vars, 0).toMat3d().sym();

	// collect all terms
	double D[6][6] = { 0 };
	D[0][0] = 2.0*dSdC00.xx(); D[0][1] = 2.0*dSdC11.xx(); D[0][2] = 2.0*dSdC22.xx(); D[0][3] = 2.0*dSdC01.xx(); D[0][4] = 2.0*dSdC12.xx(); D[0][5] = 2.0*dSdC02.xx();
	D[1][0] = 2.0*dSdC00.yy(); D[1][1] = 2.0*dSdC11.yy(); D[1][2] = 2.0*dSdC22.yy(); D[1][3] = 2.0*dSdC01.yy(); D[1][4] = 2.0*dSdC12.yy(); D[1][5] = 2.0*dSdC02.yy();
	D[2][0] = 2.0*dSdC00.zz(); D[2][1] = 2.0*dSdC11.zz(); D[2][2] = 2.0*dSdC22.zz(); D[2][3] = 2.0*dSdC01.zz(); D[2][4] = 2.0*dSdC12.zz(); D[2][5] = 2.0*dSdC02.zz();
	D[3][0] = 2.0*dSdC00.xy(); D[3][1] = 2.0*dSdC11.xy(); D[3][2] = 2.0*dSdC22.xy(); D[3][3] = 2.0*dSdC01.xy(); D[3][4] = 2.0*dSdC12.xy(); D[3][5] = 2.0*dSdC02.xy();
	D[4][0] = 2.0*dSdC00.yz(); D[4][1] = 2.0*dSdC11.yz(); D[4][2] = 2.0*dSdC22.yz(); D[4][3] = 2.0*dSdC01.yz(); D[4][4] = 2.0*dSdC12.yz(); D[4][5] = 2.0*dSdC02.yz();
	D[5][0] = 2.0*dSdC00.xz(); D[5][1] = 2.0*dSdC11.xz(); D[5][2] = 2.0*dSdC22.xz(); D[5][3] = 2.0*dSdC01.xz(); D[5][4] = 2.0*dSdC12.xz(); D[5][5] = 2.0*dSdC02.xz();
	tens4ds C(D);

	// push-forward to spatial frame
	tens4ds c = ep.push_forward(C);

	return c;
}

double FEScriptedElasticMaterial::StrainEnergyDensity(FEMaterialPoint& pt)
{
	assert(false);
	return 0.0;
}
