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



#include "FEThermoFluid.h"
#include "FEThermoFluidMaterialPoint.h"
#include <FEBioFluid/FEFluidMaterialPoint.h>
#include <FECore/FECoreKernel.h>
#include <FECore/DumpStream.h>
#include <FECore/FEModel.h>
#include "FELinearThermoElasticFluid.h"
// define the material parameters
BEGIN_FECORE_CLASS(FEThermoFluid, FEThermoFluidMaterial)

	// material parameters
    ADD_PARAMETER(m_k      , FE_RANGE_GREATER_OR_EQUAL(0.0), "k")->setUnits(UNIT_PRESSURE);
    ADD_PROPERTY(m_pElastic, "elastic", FEProperty::Optional)->SetLongName("Elastic response");

END_FECORE_CLASS();

//============================================================================
// FEThermoFluid
//============================================================================

//-----------------------------------------------------------------------------
//! FEFluid constructor

FEThermoFluid::FEThermoFluid(FEModel* pfem) : FEThermoFluidMaterial(pfem)
{
    m_k = 0;
    m_pElastic = nullptr;
}

//-----------------------------------------------------------------------------
//! FEFluid initialization
bool FEThermoFluid::Init()
{
    m_Tr = GetGlobalConstant("T");
    if (m_pElastic == nullptr) {
//        m_pElastic = fecore_alloc(FELinearThermoElasticFluid, GetFEModel());
        m_pElastic =  new FELinearThermoElasticFluid(GetFEModel());
    }

	// We set the parent of the elastic material to this fluid,
	// so that it can access the fluid parameters if needed.
	m_pElastic->SetParentFluid(this);

    return FEThermoFluidMaterial::Init();
}

//-----------------------------------------------------------------------------
void FEThermoFluid::Serialize(DumpStream& ar)
{
    FEThermoFluidMaterial::Serialize(ar);
    if (ar.IsShallow()) return;
    
    ar & m_Tr;
}

//-----------------------------------------------------------------------------
//! returns a pointer to a new material point object
FEMaterialPointData* FEThermoFluid::CreateMaterialPointData()
{
    FEFluidMaterialPoint* fpt = new FEFluidMaterialPoint();
    FEThermoFluidMaterialPoint* tfp = new FEThermoFluidMaterialPoint(fpt);
	return tfp;
}

//-----------------------------------------------------------------------------
//! bulk modulus
double FEThermoFluid::BulkModulus(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& vt = *mp.ExtractData<FEFluidMaterialPoint>();
    return -(vt.m_ef+1)*Tangent_Pressure_Strain(mp);
}

//-----------------------------------------------------------------------------
//! elastic pressure
double FEThermoFluid::Pressure(FEMaterialPoint& mp)
{
    return m_pElastic->Pressure(mp);
}

//-----------------------------------------------------------------------------
//! elastic pressure from dilatation
double FEThermoFluid::Pressure(const double e, const double T)
{
    return m_pElastic->Pressure(e, T);
}

//-----------------------------------------------------------------------------
double FEThermoFluid::Tangent_Pressure_Strain(FEMaterialPoint& mp)
{
    return m_pElastic->Tangent_Strain(mp);
}

//-----------------------------------------------------------------------------
double FEThermoFluid::Tangent_Pressure_Strain_Strain(FEMaterialPoint& mp)
{
    return m_pElastic->Tangent_Strain_Strain(mp);
}

//-----------------------------------------------------------------------------
//! The stress of a fluid material is the sum of the fluid pressure
//! and the viscous stress.

mat3ds FEThermoFluid::Stress(FEMaterialPoint& mp)
{
	// calculate solid material stress
	mat3ds s = GetViscous()->Stress(mp);
    
    double p = Pressure(mp);
	
	// add fluid pressure
	s.xx() -= p;
	s.yy() -= p;
	s.zz() -= p;
	
	return s;
}

//-----------------------------------------------------------------------------
//! The tangent of stress with respect to strain J of a fluid material is the
//! sum of the tangent of the fluid pressure and that of the viscous stress.

mat3ds FEThermoFluid::Tangent_Strain(FEMaterialPoint& mp)
{
    // get tangent of viscous stress
    mat3ds sJ = GetViscous()->Tangent_Strain(mp);
    
    // add tangent of fluid pressure
    double dp = Tangent_Pressure_Strain(mp);
    sJ.xx() -= dp;
    sJ.yy() -= dp;
    sJ.zz() -= dp;
    
    return sJ;
}

//-----------------------------------------------------------------------------
//! calculate strain energy density (per reference volume)
double FEThermoFluid::StrainEnergyDensity(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    double sed = m_k*(fp.m_ef-log(fp.m_ef+1));
    return sed;
}

//-----------------------------------------------------------------------------
//! invert pressure-dilatation relation
bool FEThermoFluid::Dilatation(const double T, const double p, double& e)
{
    e = -p/m_k;
    return true;
}
