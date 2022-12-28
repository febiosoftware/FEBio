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
#include "FEConstPrestrain.h"

//-----------------------------------------------------------------------------
FEConstPrestrainGradient::MaterialPointData::MaterialPointData() 
{
	Fp.unit();
}
		
//-----------------------------------------------------------------------------
void FEConstPrestrainGradient::MaterialPointData::Init(bool bflag) 
{
}

//-----------------------------------------------------------------------------
FEMaterialPointData* FEConstPrestrainGradient::MaterialPointData::Copy()
{ 
	MaterialPointData* pm = new MaterialPointData(*this); 
	return pm;
}

//-----------------------------------------------------------------------------
void FEConstPrestrainGradient::MaterialPointData::Serialize(DumpStream& ar)
{
	FEMaterialPointData::Serialize(ar);
	ar & Fp;
}

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEConstPrestrainGradient, FEPrestrainGradient)
	ADD_PARAMETER(m_ramp, "ramp");
	ADD_PARAMETER(m_Fp, "F0");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEConstPrestrainGradient::FEConstPrestrainGradient(FEModel* pfem) : FEPrestrainGradient(pfem)
{
	m_ramp = 1.0;
	m_Fp = mat3dd(1.0);
}

//-----------------------------------------------------------------------------
FEMaterialPointData* FEConstPrestrainGradient::CreateMaterialPointData()
{
	return new MaterialPointData;
}

//-----------------------------------------------------------------------------
mat3d FEConstPrestrainGradient::Prestrain(FEMaterialPoint& mp)
{
	MaterialPointData& ep = *mp.ExtractData<MaterialPointData>();

	mat3d Fp = mat3dd(1.0) * (1.0 - m_ramp) + m_Fp(mp) * m_ramp;

	return Fp*ep.Fp;
}

//-----------------------------------------------------------------------------
void FEConstPrestrainGradient::Initialize(const mat3d& F, FEMaterialPoint& mp)
{
	FEElasticMaterialPoint* pt = mp.ExtractData<FEElasticMaterialPoint>();
	FEConstPrestrainGradient::MaterialPointData* ps = mp.ExtractData<FEConstPrestrainGradient::MaterialPointData>();
/*
	// calculate left polar decomposition
	mat3d R;
	mat3ds V;
	F.left_polar(V, R);

	// assign it to the pre-strain gradient
	ps->Fp = V;

	// adjust fiber vector
	vec3d a0 = pt->m_Q.col(0);
	vec3d a1 = R*a0;

	// setup orthogonal system
	vec3d b(0,0,1);
	if (b*a1 > 0.9) b = vec3d(0,1,0);
	vec3d a3 = a1^b;
	vec3d a2 = a3^a1;

	mat3d Q;
	Q(0,0) = a1.x; Q(0,1) = a2.x; Q(0,2) = a3.x;
	Q(1,0) = a1.y; Q(1,1) = a2.y; Q(1,2) = a3.y;
	Q(2,0) = a1.z; Q(2,1) = a2.z; Q(2,2) = a3.z;

	pt->m_Q = Q;
*/
	ps->Fp = F;
}
