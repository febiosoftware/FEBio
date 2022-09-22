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
#include "FEDiscreteElasticMaterial.h"

FEMaterialPointData* FEDiscreteElasticMaterialPoint::Copy()
{
	FEDiscreteElasticMaterialPoint* pt = new FEDiscreteElasticMaterialPoint(*this);
	if (m_pNext) pt->m_pNext = m_pNext->Copy();
	return pt;
}

//=============================================================================
FEDiscreteElasticMaterial::FEDiscreteElasticMaterial(FEModel* pfem) : FEDiscreteMaterial(pfem)
{

}

FEMaterialPointData* FEDiscreteElasticMaterial::CreateMaterialPointData()
{
	return new FEDiscreteElasticMaterialPoint;
}

//=============================================================================
//					FECompositeDiscreteElasticMaterial
//=============================================================================

BEGIN_FECORE_CLASS(FECompositeDiscreteMaterial, FEDiscreteElasticMaterial)
	ADD_PROPERTY(m_mats, "mat");
END_FECORE_CLASS();


FECompositeDiscreteMaterial::FECompositeDiscreteMaterial(FEModel* pfem) : FEDiscreteElasticMaterial(pfem)
{

}

vec3d FECompositeDiscreteMaterial::Force(FEDiscreteMaterialPoint& mp)
{
	vec3d force; 

	for(auto mat : m_mats)
	{
		force += mat->Force(mp);
	}

	return force;
}

mat3d FECompositeDiscreteMaterial::Stiffness(FEDiscreteMaterialPoint& mp)
{
	mat3d stiffness; stiffness.zero();

	for(auto mat : m_mats)
	{
		stiffness += mat->Stiffness(mp);
	}

	return stiffness;
}

double FECompositeDiscreteMaterial::StrainEnergy(FEDiscreteMaterialPoint& mp)
{
	double energy = 0; 

	for(auto mat : m_mats)
	{
		energy += mat->StrainEnergy(mp);
	}

	return energy;
}