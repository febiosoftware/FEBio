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

#include "FEKinematicGrowth.h"
#include <FECore/FECoreKernel.h>
#include <FECore/FEModel.h>
#include <FECore/log.h>
#include "FEUncoupledMaterial.h"

//-----------------------------------------------------------------------------
//! Material point
//
FEMaterialPointData* FEKinematicMaterialPoint::Copy()
{
    FEKinematicMaterialPoint* pt = new FEKinematicMaterialPoint(*this);
    if (m_pNext) pt->m_pNext = m_pNext->Copy();
    return pt;
}

//-----------------------------------------------------------------------------
void FEKinematicMaterialPoint::Init()
{
    FEMaterialPointData::Init();
    
    // intialize data
    m_rhor = 0;
    m_Fe = mat3dd(1);
    m_Je = 1;
    m_Fg = mat3dd(1);
}

//-----------------------------------------------------------------------------
void FEKinematicMaterialPoint::Update(const FETimeInfo& timeInfo)
{
    FEMaterialPointData::Update(timeInfo);
}

//-----------------------------------------------------------------------------
void FEKinematicMaterialPoint::Serialize(DumpStream& ar)
{
    if (ar.IsSaving())
    {
        ar << m_Fe << m_Fg << m_rhor;
    }
    else
    {
        ar >> m_Fe >> m_Fg >> m_rhor;
    }
    FEMaterialPointData::Serialize(ar);
}


//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEKinematicGrowth, FEElasticMaterial)
    ADD_PROPERTY(m_pBase, "elastic");
    ADD_PROPERTY(m_pGrowth, "growth");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FEKinematicGrowth::FEKinematicGrowth(FEModel* pfem) : FEElasticMaterial(pfem)
{
    m_pBase = nullptr;
    m_pGrowth = nullptr;
}

//-----------------------------------------------------------------------------
//! create material point
FEMaterialPointData* FEKinematicGrowth::CreateMaterialPointData()
{
    FEElasticMaterial* pme = GetBaseMaterial();
    FEMaterialPointData* ep = pme->CreateMaterialPointData();
    return new FEKinematicMaterialPoint(ep);
}

//-----------------------------------------------------------------------------
//! Get the elastic deformation
FEMaterialPoint FEKinematicGrowth::GetElasticDeformationMaterialPoint(FEMaterialPoint& mp)
{
    // Get the growth tensor inverse
    FEGrowthTensor* gmat = GetGrowthMaterial();
    // material axes
    mat3d Q = GetLocalCS(mp);
    // get the fiber vector in local coordinates
    vec3d fiber = gmat->m_fiber->unitVector(mp);
    // convert to global coordinates
    vec3d a0 = Q * fiber;
    
    mat3d Fgi = gmat->GrowthTensorInverse(mp, a0);
    double Jgi = Fgi.det();
    
    // Get the deformation gradient and evaluate elastic deformation
    FEMaterialPoint pt = mp;
    FEElasticMaterialPoint& pe = *pt.ExtractData<FEElasticMaterialPoint>();
    mat3d Fe = pe.m_F*Fgi;
    double Je = pe.m_J*Jgi;
    // substitute elastic deformation in material point
    pe.m_F = Fe;
    pe.m_J = Je;

    return pt;
}

//-----------------------------------------------------------------------------
//! data initialization
bool FEKinematicGrowth::Init()
{
    FEUncoupledMaterial* m_pMat = dynamic_cast<FEUncoupledMaterial*>((FEElasticMaterial*)m_pBase);
    if (m_pMat != nullptr) {
        feLogError("Elastic material should not be of type uncoupled");
        return false;
    }
    
    return FEElasticMaterial::Init();
}

//-----------------------------------------------------------------------------
//! Returns the Cauchy stress
mat3ds FEKinematicGrowth::Stress(FEMaterialPoint& mp)
{
    // evaluate stress
    FEElasticMaterial* emat = GetBaseMaterial();
    mat3ds s = emat->Stress(mp);

    return s;
}

//-----------------------------------------------------------------------------
//! Returns the spatial tangent
tens4ds FEKinematicGrowth::Tangent(FEMaterialPoint& mp)
{
    // evaluate tangent
    FEElasticMaterial* emat = GetBaseMaterial();
    tens4ds c = emat->Tangent(mp);

    return c;
}

//-----------------------------------------------------------------------------
//! Returns the strain energy density
double FEKinematicGrowth::StrainEnergyDensity(FEMaterialPoint& mp)
{
    // evaluate sed
    FEElasticMaterial* emat = GetBaseMaterial();
    double sed = emat->StrainEnergyDensity(mp);

    return sed;
}

//-----------------------------------------------------------------------------
//! update material point at each iteration
void FEKinematicGrowth::UpdateSpecializedMaterialPoints(FEMaterialPoint& mp, const FETimeInfo& tp)
{
    // Get the growth tensor inverse
    FEGrowthTensor* gmat = GetGrowthMaterial();
    // material axes
    mat3d Q = GetLocalCS(mp);
    // get the fiber vector in local coordinates
    vec3d fiber = gmat->m_fiber->unitVector(mp);
    // convert to global coordinates
    vec3d a0 = Q * fiber;
    
    FEElasticMaterialPoint& pe = *mp.ExtractData<FEElasticMaterialPoint>();
    
    // Get the deformation gradient and evaluate elastic deformation
    mat3d Fg = gmat->GrowthTensor(mp, a0);
    mat3d Fe = pe.m_F*gmat->GrowthTensorInverse(mp, a0);

    // extract Kinematic growth material point
    FEKinematicMaterialPoint& pt = *mp.ExtractData<FEKinematicMaterialPoint>();
    
    pt.m_Fg = Fg;
    pt.m_Fe = Fe;
    pt.m_Je = Fe.det();
    pt.m_rhor = GetBaseMaterial()->Density(mp)*GetGrowthMaterial()->GrowthDensity(mp, a0);
    
    // overwrite the elastic material point deformation gradient information
    // to match elastic deformation gradient of kinematic growth material
    pe.m_F = pt.m_Fe;
    pe.m_J = pt.m_Je;
}
