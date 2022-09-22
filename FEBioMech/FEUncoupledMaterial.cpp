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
#include "FEUncoupledMaterial.h"
#include <FECore/log.h>

//-----------------------------------------------------------------------------
// Material parameters for FEUncoupledMaterial
BEGIN_FECORE_CLASS(FEUncoupledMaterial, FEElasticMaterial)
	ADD_PARAMETER(m_K      , FE_RANGE_GREATER_OR_EQUAL(0.0), "k")->setUnits(UNIT_PRESSURE)->MakeTopLevel(true)->setLongName("bulk modulus");
    ADD_PARAMETER(m_npmodel, "pressure_model")->setEnums("default\0NIKE3D\0Abaqus\0Abaqus (GOH)\0")->MakeTopLevel(true);
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FEUncoupledMaterial::FEUncoupledMaterial(FEModel* pfem) : FEElasticMaterial(pfem)
{
	m_K = 0;	// invalid value!
    m_npmodel = 0;
}

//-----------------------------------------------------------------------------
bool FEUncoupledMaterial::Init()
{
	// K has to be non-zero if this is a top-level material, otherwise
	// it should be zero. 
	FECoreBase* parent = GetParent();
	FEUncoupledMaterial* parentMat = dynamic_cast<FEUncoupledMaterial*>(GetParent());
	if (parentMat)
	{
		if (m_K != 0.0)
		{
			feLogWarning("K should only be defined at the top-level.\nValue will be ignored.");

			// NOTE: This should not be neccessary, but just to be safe. 
			m_K = 0.0;
		}
	}
	else
	{
		if (m_K <= 0.0)
		{
			feLogError("K must be a positive number.");
			return false;
		}
	}
	return FEElasticMaterial::Init();
}

//-----------------------------------------------------------------------------
//! The stress function calculates the total Cauchy stress as a sum of 
//! two terms, namely the deviatoric stress and the pressure. 
mat3ds FEUncoupledMaterial::Stress(FEMaterialPoint &mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// calculate the stress as a sum of deviatoric stress and pressure
    pt.m_p = UJ(pt.m_J);
	return mat3dd(pt.m_p) + DevStress(mp);
}

//------------------------------------------------------------------------------
//! The tangent function calculates the total spatial tangent, that is it calculates
//! the push-forward of the derivative of the 2ndPK stress with respect to C. However,
//! for an uncoupled material, the 2ndPK stress decouples in a deviatoric and a 
//! dilatational component. The deviatoric tangent is provided by the particular
//! material and the dilatational component is added here.
//!
tens4ds FEUncoupledMaterial::Tangent(FEMaterialPoint &mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// 2nd-order identity tensor
	mat3dd I(1);

	// 4th-order identity tensors
	tens4ds IxI = dyad1s(I);
	tens4ds I4  = dyad4s(I);
	
	// pressure
	pt.m_p = UJ(pt.m_J);
	
	// tangent is sum of three terms
	// C = c_tilde + c_pressure + c_k
	//
	// + c_tilde is the derivative of the deviatoric stress with respect to C
	// + c_pressure is p*d(JC)/dC
	// + c_k comes from the derivative of p with respect to C
	// 
	// Note that the c_k term is not necessary in the 3F formulation (since p is independant variable) 
	// but we do need to add it here.
	//
	//        c_tilde         c_pressure            c_k
	return DevTangent(mp) + (IxI - I4*2)*pt.m_p + IxI*(UJJ(pt.m_J)*pt.m_J);
}

//-----------------------------------------------------------------------------
//! The strain energy density function calculates the total sed as a sum of
//! two terms, namely the deviatoric sed and U(J).
double FEUncoupledMaterial::StrainEnergyDensity(FEMaterialPoint &mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
	// calculate the stress as a sum of deviatoric stress and pressure
	return U(pt.m_J) + DevStrainEnergyDensity(mp);
}

//-----------------------------------------------------------------------------
//! The strain energy density function calculates the total sed as a sum of
//! two terms, namely the deviatoric sed and U(J).
double FEUncoupledMaterial::StrongBondSED(FEMaterialPoint &mp)
{
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // calculate the stress as a sum of deviatoric stress and pressure
    return U(pt.m_J) + StrongBondDevSED(mp);
}

//-----------------------------------------------------------------------------
//! The strain energy density function calculates the total sed as a sum of
//! two terms, namely the deviatoric sed and U(J).
double FEUncoupledMaterial::WeakBondSED(FEMaterialPoint &mp)
{
    // calculate the stress as a sum of deviatoric stress and pressure
    return WeakBondDevSED(mp);
}

//-----------------------------------------------------------------------------
FEMaterialPointData* FEUncoupledMaterial::CreateMaterialPointData()
{
	FEElasticMaterialPoint* mp = new FEElasticMaterialPoint;
	mp->m_buncoupled = true;
	return mp;
}
