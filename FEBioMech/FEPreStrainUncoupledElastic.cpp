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
#include "FEPreStrainUncoupledElastic.h"

BEGIN_FECORE_CLASS(FEPreStrainUncoupledElastic, FEUncoupledMaterial)
	ADD_PROPERTY(m_mat, "elastic");
	ADD_PROPERTY(m_Fp, "prestrain", false);
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEPreStrainUncoupledElastic::FEPreStrainUncoupledElastic(FEModel* pfem) : FEUncoupledMaterial(pfem)
{
	m_mat = nullptr;
	m_Fp = nullptr;
}

//-----------------------------------------------------------------------------
//! Create material point data for this material
FEMaterialPointData* FEPreStrainUncoupledElastic::CreateMaterialPointData()
{ 
	FEMaterialPointData* pm = m_mat->CreateMaterialPointData();
	if (m_Fp) pm->Append(m_Fp->CreateMaterialPointData());
	return new FEPrestrainMaterialPoint(pm);
}

//-----------------------------------------------------------------------------
//! calculate (pre-strained) density 
double FEPreStrainUncoupledElastic::Density(FEMaterialPoint& mp)
{
	double d0 = FEElasticMaterial::Density(mp);

	mat3d Fp = PrestrainGradient(mp);
	double Jp = Fp.det();

	return d0 / Jp;
}

//-----------------------------------------------------------------------------
mat3d FEPreStrainUncoupledElastic::PrestrainGradient(FEMaterialPoint& mp)
{
	mat3d F0 = mat3d::identity();
	if (m_Fp) F0 = m_Fp->Prestrain(mp);

	FEPrestrainMaterialPoint& pt = *(mp.ExtractData<FEPrestrainMaterialPoint>());
	pt.setInitialPrestrain(F0);

	return pt.prestrain();
}

//-----------------------------------------------------------------------------
mat3ds FEPreStrainUncoupledElastic::DevStress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& ep = *(mp.ExtractData<FEElasticMaterialPoint>());
	FEPrestrainMaterialPoint& pp = *(mp.ExtractData<FEPrestrainMaterialPoint>());

	// store the original deformation gradient
	mat3d F0 = ep.m_F;
	double J0 = ep.m_J;

	// get the pre-strain deformation gradient
	mat3d Fp = PrestrainGradient(mp);

	// pre-multiply the pre-strain
	ep.m_F = ep.m_F*Fp;
	ep.m_J = ep.m_F.det();

	// evaluate the stress
	mat3ds s = m_mat->DevStress(mp);

	// restore original deformation gradient
	ep.m_F = F0;
	ep.m_J = J0;

	// return stress
	return s;
}

//-----------------------------------------------------------------------------
tens4ds FEPreStrainUncoupledElastic::DevTangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& ep = *(mp.ExtractData<FEElasticMaterialPoint>());
	FEPrestrainMaterialPoint& pt = *(mp.ExtractData<FEPrestrainMaterialPoint>());

	// store the original deformation gradient
	mat3d F0 = ep.m_F;
	double J0 = ep.m_J;

	// get the pre-strain deformation gradient
	mat3d Fp = pt.prestrain();

	// pre-multiply the pre-strain
	ep.m_F = ep.m_F*Fp;
	ep.m_J = ep.m_F.det();

	// evaluate the tangent
	tens4ds c = m_mat->DevTangent(mp);

	// restore original deformation gradient
	ep.m_F = F0;
	ep.m_J = J0;

	// return spatial tangent
	return c;
}
