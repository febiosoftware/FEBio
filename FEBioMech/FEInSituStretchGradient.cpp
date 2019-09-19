/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in
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
#include "FEInSituStretchGradient.h"

//-----------------------------------------------------------------------------
FEInSituStretchGradient::MaterialPointData::MaterialPointData()
{
	m_lam = 0.0;
}
	
void FEInSituStretchGradient::MaterialPointData::Init(bool bflag)
{
}

FEMaterialPoint* FEInSituStretchGradient::MaterialPointData::Copy()
{
	MaterialPointData* pm = new MaterialPointData();
	pm->m_lam = m_lam;
	return pm;
}

void FEInSituStretchGradient::MaterialPointData::Serialize(DumpStream& ar)
{
	if (ar.IsSaving())
	{
		ar << m_lam;
	}
	else
	{
		ar >> m_lam;
	}
}

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEInSituStretchGradient, FEPrestrainGradient)
	ADD_PARAMETER(m_lam , "stretch"  );
	ADD_PARAMETER(m_biso, "isochoric");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEInSituStretchGradient::FEInSituStretchGradient(FEModel* pfem) : FEPrestrainGradient(pfem) 
{
	m_lam = 1.0;
	m_biso = true;
}

//-----------------------------------------------------------------------------
FEMaterialPoint* FEInSituStretchGradient::CreateMaterialPointData()
{
	return new MaterialPointData;
}

//-----------------------------------------------------------------------------
double FEInSituStretchGradient::InSituStretch(FEMaterialPoint& mp)
{
	MaterialPointData& psp = *mp.ExtractData<MaterialPointData>();
	// get the target stretch
	double trg = 0;
	if (psp.m_lam ==0.0) trg = (m_lam < 1e-15 ? 1.0 : m_lam);
	else
	{
		double w = m_lam;
		double l = psp.m_lam;
		trg = w*l + (1.0-w);
	}
	return trg;
}

//-----------------------------------------------------------------------------
mat3d FEInSituStretchGradient::Prestrain(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	MaterialPointData& psp = *mp.ExtractData<MaterialPointData>();

	// get the in-situ stretch
	double lam = InSituStretch(mp);

	// set-up local uni-axial stretch tensor
	double l = lam;
	double li = (m_biso ? 1.0/sqrt(l) : 1.0);
	mat3d U(l, 0.0, 0.0, 0.0, li, 0.0, 0.0, 0.0, li);

	// get the coordinate transformation
	mat3d Q = GetLocalCS(mp);
	mat3d Qt = Q.transpose();

	// rotate stretch tensor to global coordinates
	mat3d F_bar = Q*U*Qt;

	return F_bar;
}

//-----------------------------------------------------------------------------
void FEInSituStretchGradient::Initialize(const mat3d& F, FEMaterialPoint& mp)
{
	FEElasticMaterialPoint* pt = mp.ExtractData<FEElasticMaterialPoint>();
	FEInSituStretchGradient::MaterialPointData* ptis = mp.ExtractData<FEInSituStretchGradient::MaterialPointData>();

	// calculate left polar decomposition
	mat3d R;
	mat3ds V;
	F.left_polar(V, R);

	// get the fiber vector
	mat3d Q = GetLocalCS(mp);
	vec3d a0 = Q*vec3d(1,0,0);
	vec3d a1 = F*a0;

	// calculate the fiber stretch
	double lam = a1.unit();

	// assign it to the pre_stretch
	ptis->m_lam = lam;

	// setup orthogonal system
	vec3d b(0,0,1);
	if (b*a1 > 0.9) b = vec3d(0,1,0);
	vec3d a3 = a1^b;
	vec3d a2 = a3^a1;

	Q(0,0) = a1.x; Q(0,1) = a2.x; Q(0,2) = a3.x;
	Q(1,0) = a1.y; Q(1,1) = a2.y; Q(1,2) = a3.y;
	Q(2,0) = a1.z; Q(2,1) = a2.z; Q(2,2) = a3.z;
//	pt->m_Q = Q;
}
