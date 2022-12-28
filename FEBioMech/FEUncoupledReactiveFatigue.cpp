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
#include "FEUncoupledReactiveFatigue.h"
#include "FEDamageCriterion.h"
#include "FEDamageCDF.h"
#include <FECore/FECoreKernel.h>
#include <FECore/DumpStream.h>
#include <FECore/log.h>
#include <FECore/FEAnalysis.h>
#include <FECore/FEModel.h>

//////////////////////////// FATIGUE MATERIAL /////////////////////////////////
//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEUncoupledReactiveFatigue, FEUncoupledMaterial)
ADD_PARAMETER(m_k0   , FE_RANGE_GREATER_OR_EQUAL(0.0), "k0"  );
ADD_PARAMETER(m_beta , FE_RANGE_GREATER_OR_EQUAL(0.0), "beta");

// set material properties
ADD_PROPERTY(m_pBase, "elastic");
ADD_PROPERTY(m_pIdmg, "elastic_damage");
ADD_PROPERTY(m_pFdmg, "fatigue_damage");
ADD_PROPERTY(m_pIcrt, "elastic_criterion");
ADD_PROPERTY(m_pFcrt, "fatigue_criterion");

END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor.
FEUncoupledReactiveFatigue::FEUncoupledReactiveFatigue(FEModel* pfem) : FEUncoupledMaterial(pfem)
{
    m_pBase = 0;
    m_pIdmg = 0;
    m_pFdmg = 0;
    m_pIcrt = 0;
    m_pFcrt = 0;
}

//-----------------------------------------------------------------------------
FEMaterialPointData* FEUncoupledReactiveFatigue::CreateMaterialPointData()
{
	return new FEReactiveFatigueMaterialPoint(m_pBase->CreateMaterialPointData());
}

//-----------------------------------------------------------------------------
//! Initialization.
bool FEUncoupledReactiveFatigue::Init()
{
    return FEUncoupledMaterial::Init();
}

//-----------------------------------------------------------------------------
//! calculate stress at material point
mat3ds FEUncoupledReactiveFatigue::DevStress(FEMaterialPoint& pt)
{
    // evaluate the damage
    double d = Damage(pt);
    
    // evaluate the stress
    mat3ds s = m_pBase->DevStress(pt);
    
    // return damaged stress
    return s*(1-d);
}

//-----------------------------------------------------------------------------
//! calculate tangent stiffness at material point
tens4ds FEUncoupledReactiveFatigue::DevTangent(FEMaterialPoint& pt)
{
    // evaluate the damage
    double d = Damage(pt);
    
    // evaluate the tangent
    tens4ds c = m_pBase->DevTangent(pt);
    
    // return damaged tangent
    return c*(1-d);
}

//-----------------------------------------------------------------------------
//! calculate strain energy density at material point
double FEUncoupledReactiveFatigue::DevStrainEnergyDensity(FEMaterialPoint& pt)
{
    // evaluate the damage
    double d = Damage(pt);
    
    // evaluate the strain energy density
    double sed = m_pBase->DevStrainEnergyDensity(pt);
    
    // return damaged sed
    return sed*(1-d);
}

//-----------------------------------------------------------------------------
//! calculate damage at material point
double FEUncoupledReactiveFatigue::Damage(FEMaterialPoint& pt)
{
    // get the reactive fatigue material point data
    FEReactiveFatigueMaterialPoint& pd = *pt.ExtractData<FEReactiveFatigueMaterialPoint>();
    
    return pd.m_D;
}

//-----------------------------------------------------------------------------
// update fatigue material point at each iteration
void FEUncoupledReactiveFatigue::UpdateSpecializedMaterialPoints(FEMaterialPoint& pt, const FETimeInfo& tp)
{
    double dt = tp.timeIncrement;
    
    // get the fatigue material point data
    FEReactiveFatigueMaterialPoint& pd = *pt.ExtractData<FEReactiveFatigueMaterialPoint>();
    
    // get damage criterion for intact bonds at current time
    pd.m_Xitrl = m_pIcrt->DamageCriterion(pt);
    if (pd.m_Xitrl > pd.m_Ximax)
        pd.m_Fit = m_pIdmg->cdf(pt,pd.m_Xitrl);
    else
        pd.m_Fit = pd.m_Fip;
    
    // get damage criterion for fatigue bonds at current time
    double Xftrl = m_pFcrt->DamageCriterion(pt);
    for (int ig=0; ig < pd.m_fb.size(); ++ig) {
        if (Xftrl > pd.m_fb[ig].m_Xfmax) {
            pd.m_fb[ig].m_Xftrl = Xftrl;
            pd.m_fb[ig].m_Fft = m_pFdmg->cdf(pt,pd.m_fb[ig].m_Xftrl);
        }
        else
            pd.m_fb[ig].m_Fft = pd.m_fb[ig].m_Ffp;
    }
    
    // evaluate time derivative of intact bond criterion
    pd.m_aXit = (pd.m_Xitrl - pd.m_Xip)/dt;
    
    // solve for bond mass fractions iteratively
    double eps = 1e-6;
    int maxit = 10;
    double wbi = 0;
    int iter = 0;
    do {
        wbi = pd.m_wbt;
        // evaluate mass supply from fatigue of intact bonds
        double dwf = m_k0*dt/2*(pow(fabs(pd.m_aXit)*pd.m_wbt,m_beta)*pd.m_wit+pow(fabs(pd.m_aXip)*pd.m_wbp,m_beta)*pd.m_wip);
        double Fdwf = m_pFdmg->cdf(pt,Xftrl);
        
        // kinetics of intact bonds
        pd.m_wit = (pd.m_Fip < 1) ? pd.m_wip*(1-pd.m_Fit)/(1-pd.m_Fip) - dwf : 0;
        pd.m_wbt = (pd.m_Fip < 1) ? pd.m_wbp + pd.m_wip*(pd.m_Fit - pd.m_Fip)/(1-pd.m_Fip) : pd.m_wbp;
        // add or update new generation
        if ((pd.m_fb.size() == 0) || pd.m_fb.back().m_time < tp.currentTime) {
            // add generation of fatigued bonds
            FatigueBond fb;
            fb.m_Fft = Fdwf;
            fb.m_wfp = dwf;
            fb.m_Xftrl = Xftrl;
            fb.m_time = tp.currentTime;
            pd.m_fb.push_back(fb);
        }
        else {
            pd.m_fb.back().m_Fft = Fdwf;
            pd.m_fb.back().m_wfp = dwf;
            pd.m_fb.back().m_Xftrl = Xftrl;
        }
        // damage kinetics of fatigued bonds
        for (int ig=0; ig < pd.m_fb.size(); ++ig) {
            pd.m_fb[ig].m_wft = (pd.m_fb[ig].m_Ffp < 1) ? pd.m_fb[ig].m_wfp*(1-pd.m_fb[ig].m_Fft)/(1-pd.m_fb[ig].m_Ffp) : 0;
            pd.m_wbt += (pd.m_fb[ig].m_Ffp < 1) ? pd.m_fb[ig].m_wfp*(pd.m_fb[ig].m_Fft-pd.m_fb[ig].m_Ffp)/(1-pd.m_fb[ig].m_Ffp) : pd.m_fb[ig].m_wfp;
        }
        // roundoff corrections
        if (pd.m_wit < 0) pd.m_wit = 0;
        if (pd.m_wbt > 1) pd.m_wbt = 1;
        // evaluate fatigue bond fraction
        pd.m_wft = 0;
        for (int ig=0; ig < pd.m_fb.size(); ++ig) pd.m_wft += pd.m_fb[ig].m_wft;
    } while ((pd.m_wbt > 0) && (pd.m_wbt < 1) && (fabs(wbi-pd.m_wbt) > eps*pd.m_wbt) && (++iter<maxit));
    
    pd.m_D = pd.m_wbt;
}
