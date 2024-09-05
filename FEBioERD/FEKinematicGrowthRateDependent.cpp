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

#include "FEKinematicGrowthRateDependent.h"
#include <FECore/FECoreKernel.h>
#include <FECore/FEModel.h>
#include <FECore/log.h>
#include <FEBioMech/FEUncoupledMaterial.h>
#include "FEGrowthTensorERD.h"
#include "FEElasticReactionDiffusion.h"
#include <FECore/tools.h>

//============================================================================
//FEKinematicMaterialPoint::FEKinematicMaterialPoint(FEMaterialPointData* ppt) : FEMaterialPointData(ppt) {}

//-----------------------------------------------------------------------------
//! Material point
//
FEMaterialPointData* FEKinematicMaterialPointERD::Copy()
{
    FEKinematicMaterialPointERD* pt = new FEKinematicMaterialPointERD(*this);
    if (m_pNext) pt->m_pNext = m_pNext->Copy();
    return pt;
}

//-----------------------------------------------------------------------------
void FEKinematicMaterialPointERD::Init()
{   
    // intialize data
    m_Fe = mat3dd(1.0);
    m_Fg = mat3dd(1.0);
    m_Je = 1.0;
    m_Jg = 1.0;
    m_theta = 1.0;
    m_theta_p = 1.0;
    m_rhor = 0;
    FEMaterialPointData::Init();
}

//-----------------------------------------------------------------------------
void FEKinematicMaterialPointERD::Update(const FETimeInfo& timeInfo)
{
    m_theta_p = m_theta;
    FEMaterialPointData::Update(timeInfo);
}

//-----------------------------------------------------------------------------
void FEKinematicMaterialPointERD::Serialize(DumpStream& ar)
{
    if (ar.IsSaving())
    {
        ar << m_Fe << m_Fg << m_Je << m_Jg << m_theta << m_theta_p << m_rhor;
    }
    else
    {
        ar >> m_Fe >> m_Fg >> m_Je >> m_Jg >> m_theta << m_theta_p >> m_rhor;
    }
    FEMaterialPointData::Serialize(ar);
}


//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEKinematicGrowthRateDependent, FEElasticMaterial)
    ADD_PROPERTY(m_pBase, "elastic");
    ADD_PROPERTY(m_pGrowth, "growth");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FEKinematicGrowthRateDependent::FEKinematicGrowthRateDependent(FEModel* pfem) : FEElasticMaterial(pfem)
{
    m_pBase = nullptr;
    m_pGrowth = nullptr;
}

//-----------------------------------------------------------------------------
//! create material point
FEMaterialPointData* FEKinematicGrowthRateDependent::CreateMaterialPointData()
{
    FEGrowthTensorERD* pmf = GetGrowthMaterial();
    FEMaterialPointData* mp = pmf->CreateMaterialPointData();
    FEMaterialPointData* ep = new FEElasticMaterialPoint();
    ep->SetNext(mp);
    return new FEKinematicMaterialPointERD(ep);
}

//-----------------------------------------------------------------------------
//! data initialization
bool FEKinematicGrowthRateDependent::Init()
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
mat3ds FEKinematicGrowthRateDependent::Stress(FEMaterialPoint& mp)
{
    // Get the growth tensor inverse
    FEGrowthTensorERD* gmat = GetGrowthMaterial();
    // material axes
    mat3d Q = GetLocalCS(mp);
    // get the fiber vector in local coordinates
    vec3d fiber = gmat->m_fiber.unitVector(mp);
    // convert to global coordinates
    vec3d a0 = Q * fiber;
    
    mat3d Fgi = gmat->GrowthTensorInverse(mp, a0);
    double Jgi = Fgi.det();
    // Get the deformation gradient and evaluate elastic deformation
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    mat3d Fe = pt.m_F * Fgi;
    double Je = pt.m_J * Jgi;
    // keep safe copy of deformation gradient.
    mat3d F = pt.m_F;
    double J = pt.m_J;
    // substitute elastic deformation in material point
    pt.m_F = Fe;
    pt.m_J = Je;
    pt.m_J = Je;
    // evaluate stress
    FEElasticMaterial* emat = GetBaseMaterial();
    // The base class stress function divides by Je rather than J when pushing forward. Thus the cauchy stress is the base class stress divided by Jg.
    mat3ds s = emat->Stress(mp);
    // restore safe copy
    pt.m_F = F;
    pt.m_J = J;

    return s;
}

//-----------------------------------------------------------------------------
//! Returns the spatial elastic tangent
tens4ds FEKinematicGrowthRateDependent::Tangent(FEMaterialPoint& mp)
{
    // Get the growth tensor inverse
    FEGrowthTensorERD* gmat = GetGrowthMaterial();
    // material axes
    mat3d Q = GetLocalCS(mp);
    // get the fiber vector in local coordinates
    vec3d fiber = gmat->m_fiber.unitVector(mp);
    // convert to global coordinates
    vec3d a0 = Q * fiber;
    
    mat3d Fgi = gmat->GrowthTensorInverse(mp, a0);
    double Jgi = Fgi.det();

    // Get the deformation gradient and evaluate elastic deformation
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    mat3d Fe = pt.m_F * Fgi;
    double Je = pt.m_J * Jgi;

    mat3d Fei = Fe.inverse();
    mat3d FeiT = Fei.transpose();

    // keep safe copy of deformation gradient.
    mat3d F = pt.m_F;
    double J = pt.m_J;
    // substitute elastic deformation in material point
    pt.m_F = Fe;
    pt.m_J = Je;
    FEElasticMaterial* emat = GetBaseMaterial();
    // ce_ijkl = (1/Je) Le_MNPQ Fe_iM Fe_jN Fe_kP Fe_lQ
    tens4ds ce = emat->Tangent(mp);
    // restore deformation gradient in material point
    pt.m_F = F;
    pt.m_J = J;

    return ce;
}

//-----------------------------------------------------------------------------
//! Returns the spatial tangent
mat3ds FEKinematicGrowthRateDependent::dSdtheta(FEMaterialPoint& mp)
{   
    //Get common variables
    // Get the growth tensor inverse
    FEGrowthTensorERD* gmat = GetGrowthMaterial();
    // material axes
    mat3d Q = GetLocalCS(mp);
    // get the fiber vector in local coordinates
    vec3d fiber = gmat->m_fiber.unitVector(mp);
    // convert to global coordinates
    vec3d a0 = Q * fiber;
    
    return dSdFg(mp).dot(gmat->dFgdtheta(mp, a0));
}

//-----------------------------------------------------------------------------
//! Returns the spatial tangent
tens4ds FEKinematicGrowthRateDependent::dSdFg(FEMaterialPoint& mp)
{
    //Get common variables
    // Get the growth tensor inverse
    FEGrowthTensorERD* gmat = GetGrowthMaterial();
    // material axes
    mat3d Q = GetLocalCS(mp);
    // get the fiber vector in local coordinates
    vec3d fiber = gmat->m_fiber.unitVector(mp);
    // convert to global coordinates
    vec3d a0 = Q * fiber;

    mat3d Fg = gmat->GrowthTensor(mp, a0);
    double Jg = Fg.det();
    mat3d Fgi = gmat->GrowthTensorInverse(mp, a0);

    // Get the deformation gradient and evaluate elastic deformation
    FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();
    FEKinematicMaterialPointERD* kp = ep.ExtractData<FEKinematicMaterialPointERD>();
    double theta = std::max(kp->m_theta, 1.0e-20);

    mat3d Fe = ep.m_F * Fgi;
    mat3d Fei = Fe.inverse();
    mat3d FeiT = Fei.transpose();
    double Je = Fe.det();
    mat3ds Ce = (Fe.transpose() * Fe).sym();
    double J = ep.m_J;
    FEElasticMaterial* emat = GetBaseMaterial();
    mat3d Fi = ep.m_F.inverse();
    // The base class stress function divides by Je rather than J when pushing forward. Thus the cauchy stress is the base class stress divided by Jg.
    mat3ds S = J * (Fi * emat->Stress(mp) * Fi.transpose()).sym();
    double dt = this->CurrentTimeIncrement();

    // solve dS/dFg
    tens4ds Ce4d = Tangent(mp);
    tens4ds Le = Je * Ce4d.pp(FeiT);
    tens4d SoFgiT = dyad1(S, Fgi.transpose());
    tens4d FgioslashS = dyad2(Fgi, S);
    tens4d SobslashFgi = dyad3(S, Fgi);
    tens4d FgioslashFgi = dyad2(Fgi, Fgi);
    tens4d CeoslashFgiT = dyad2(Ce, Fgi.transpose());
    tens4d FgiTobslashCe = dyad3(Fgi.transpose(), Ce);

    tens4ds dSdFg = (SoFgiT -FgioslashS - SobslashFgi - 0.5 * Jg * ddot(FgioslashFgi, ddot(Le, (FgiTobslashCe + CeoslashFgiT)))).supersymm();
    return dSdFg;
}

//-----------------------------------------------------------------------------
//! Returns the spatial tangent
mat3ds FEKinematicGrowthRateDependent::dTdc(FEMaterialPoint& mp, int sol)
{
    mat3ds dTdc = mat3ds(0.0);
    FEGrowthTensorERD* gmat = GetGrowthMaterial();
    if (sol == gmat->m_sol_id - 1)
    {
        // Get the deformation gradient and evaluate elastic deformation
        FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();
        mat3d F = ep.m_F;
        double J = ep.m_J;
        double k = gmat->ActivationFunction(mp);
        double g = gmat->m_gm(mp);
        mat3d Q = GetLocalCS(mp);
        // get the fiber vector in local coordinates
        vec3d fiber = gmat->m_fiber.unitVector(mp);
        // convert to global coordinates
        vec3d a0 = Q * fiber;

        dTdc = (F * dSdtheta(mp) * F.transpose()).sym() * (k * g / J);
    }
    return dTdc;
}

//-----------------------------------------------------------------------------
//! Returns the strain energy density
double FEKinematicGrowthRateDependent::StrainEnergyDensity(FEMaterialPoint& mp)
{
    // Get the growth tensor inverse
    FEGrowthTensorERD* gmat = GetGrowthMaterial();
    // material axes
    mat3d Q = GetLocalCS(mp);
    // get the fiber vector in local coordinates
    vec3d fiber = gmat->m_fiber.unitVector(mp);
    // convert to global coordinates
    vec3d a0 = Q * fiber;

    mat3d Fgi = gmat->GrowthTensorInverse(mp, a0);
    double Jgi = Fgi.det();

    // Get the deformation gradient and evaluate elastic deformation
    FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    mat3d Fe = pt.m_F * Fgi;
    double Je = pt.m_J * Jgi;
    // keep safe copy of deformation gradient.
    mat3d F = pt.m_F;
    double J = pt.m_J;
    // substitute elastic deformation in material point
    pt.m_F = Fe;
    pt.m_J = Je;
    // evaluate stress
    FEElasticMaterial* emat = GetBaseMaterial();
    double sed = emat->StrainEnergyDensity(mp);
    // restore safe copy
    pt.m_F = F;
    pt.m_J = J;

    return sed;
}

//-----------------------------------------------------------------------------
//! update material point at each iteration
void FEKinematicGrowthRateDependent::UpdateSpecializedMaterialPoints(FEMaterialPoint& mp, const FETimeInfo& tp)
{
    // Get the growth tensor inverse
    FEGrowthTensorERD* gmat = GetGrowthMaterial();
    // material axes
    mat3d Q = GetLocalCS(mp);
    // get the fiber vector in local coordinates
    vec3d fiber = gmat->m_fiber.unitVector(mp);
    // convert to global coordinates
    vec3d a0 = Q * fiber;

    // extract Elastic growth material point
    FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();
    // extract Kinematic growth material point
    FEKinematicMaterialPointERD& kp = *mp.ExtractData<FEKinematicMaterialPointERD>();
    kp.m_Fg = gmat->GrowthTensor(mp, a0);
    kp.m_Fe = ep.m_F * kp.m_Fg.inverse();
    kp.m_Je = kp.m_Fe.det();
    kp.m_Jg = kp.m_Fg.det();
    double dt = GetTimeInfo().timeIncrement;
    kp.m_theta = kp.m_theta_p + gmat->GrowthRate(mp) * dt;
    kp.m_rhor = GetBaseMaterial()->Density(mp)*GetGrowthMaterial()->GrowthDensity(mp, a0);
}
