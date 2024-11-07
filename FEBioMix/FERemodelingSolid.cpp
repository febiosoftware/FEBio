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
#include "FERemodelingSolid.h"
#include <FECore/log.h>
#include <FEBioMech/FEElasticMixture.h>

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FERemodelingSolid, FEElasticMaterial)
	ADD_PARAMETER(m_sbm , "sbm")->setEnums("$(sbms)");
    ADD_PROPERTY (m_pMat, "solid");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
bool FERemodelingSolid::Init()
{
	if (FEElasticMaterial::Init() == false) return false;
	
	// get the parent material which must be a multiphasic material
    m_pMP = GetAncestor()->ExtractProperty<FEMultiphasic>();
	if (m_pMP == 0) {
		feLogError("Parent material must be multiphasic");
		return false;
	}

	// extract the local id of the SBM whose density controls Young's modulus from the global id
	m_lsbm = m_pMP->FindLocalSBMID(m_sbm);
	if (m_lsbm == -1) {
		feLogError("Invalid value for sbm");
		return false;
	}

    FEElasticMaterial* pem = m_pMP->GetSolid();
    FEElasticMixture* psm = dynamic_cast<FEElasticMixture*>(pem);
    if (psm == nullptr) {
        m_comp = -1;    // in case material is not a solid mixture
        return true;
    }
    
    for (int i=0; i<psm->Materials(); ++i) {
        pem = psm->GetMaterial(i);
        if (pem == this) {
            m_comp = i;
            break;
        }
    }

	return true;
}

//-----------------------------------------------------------------------------
void FERemodelingSolid::Serialize(DumpStream& ar)
{
	FEElasticMaterial::Serialize(ar);
	if (ar.IsShallow()) return;
	ar & m_lsbm;
}

//-----------------------------------------------------------------------------
//! Create material point data
FEMaterialPointData* FERemodelingSolid::CreateMaterialPointData()
{
	return new FERemodelingMaterialPoint(new FEElasticMaterialPoint);
}

//-----------------------------------------------------------------------------
//! update specialize material point data
void FERemodelingSolid::UpdateSpecializedMaterialPoints(FEMaterialPoint& mp, const FETimeInfo& tp)
{
    FESolutesMaterialPoint& spt = *mp.ExtractData<FESolutesMaterialPoint>();
    FERemodelingMaterialPoint* rpt = mp.ExtractData<FERemodelingMaterialPoint>();
    rpt->m_rhor = spt.m_sbmr[m_lsbm];
    rpt->m_rhorp = spt.m_sbmrp[m_lsbm];
    rpt->m_sed = StrainEnergyDensity(mp);
}

//-----------------------------------------------------------------------------
double FERemodelingSolid::StrainEnergy(FEMaterialPoint& mp)
{
	FESolutesMaterialPoint& spt = *mp.ExtractData<FESolutesMaterialPoint>();
	
    double rhor = spt.m_sbmr[m_lsbm];
    double rho0 = m_pMP->SBMDensity(m_lsbm);
    double w = rhor/rho0;
    double sed = m_pMat->StrainEnergyDensity(mp)*w;
	
	return sed;
}

//-----------------------------------------------------------------------------
mat3ds FERemodelingSolid::Stress(FEMaterialPoint& mp)
{
	FESolutesMaterialPoint& spt = *mp.ExtractData<FESolutesMaterialPoint>();
	
    double rhor = spt.m_sbmr[m_lsbm];
    double rho0 = m_pMP->SBMDensity(m_lsbm);
    double w = rhor/rho0;
    mat3ds s = m_pMat->Stress(mp)*w;
	
	return s;
}

//-----------------------------------------------------------------------------
tens4ds FERemodelingSolid::Tangent(FEMaterialPoint& mp)
{
	FESolutesMaterialPoint& spt = *mp.ExtractData<FESolutesMaterialPoint>();
	
    double rhor = spt.m_sbmr[m_lsbm];
    double rho0 = m_pMP->SBMDensity(m_lsbm);
    double w = rhor/rho0;
    tens4ds c = m_pMat->Tangent(mp)*w;
	return c;
}

//-----------------------------------------------------------------------------
//! evaluate referential mass density
double FERemodelingSolid::Density(FEMaterialPoint& pt)
{
    FERemodelingMaterialPoint* rpt = pt.ExtractData<FERemodelingMaterialPoint>();
    if (rpt) return rpt->m_rhor;
    else {
        FEElasticMixtureMaterialPoint* emp = pt.ExtractData<FEElasticMixtureMaterialPoint>();
        if (emp) {
            rpt = emp->GetPointData(m_comp)->ExtractData<FERemodelingMaterialPoint>();
            if (rpt) return rpt->m_rhor;
        }
    }
    return 0.0;
}

//-----------------------------------------------------------------------------
//! calculate tangent of strain energy density with mass density
double FERemodelingSolid::Tangent_SE_Density(FEMaterialPoint& mp)
{
    double rho0 = m_pMP->SBMDensity(m_lsbm);
    return m_pMat->StrainEnergyDensity(mp)/rho0;
}

//-----------------------------------------------------------------------------
//! calculate tangent of stress with mass density
mat3ds FERemodelingSolid::Tangent_Stress_Density(FEMaterialPoint& mp)
{
    double rho0 = m_pMP->SBMDensity(m_lsbm);
    return m_pMat->Stress(mp)/rho0;
}

