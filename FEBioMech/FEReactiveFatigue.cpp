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
#include "FEReactiveFatigue.h"
#include "FEDamageCriterion.h"
#include "FEDamageCDF.h"
#include "FEUncoupledMaterial.h"
#include <FECore/FECoreKernel.h>
#include <FECore/DumpStream.h>
#include <FECore/log.h>
#include <FECore/FEAnalysis.h>
#include <FECore/FEModel.h>

//////////////////////////// FATIGUE MATERIAL /////////////////////////////////
//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEReactiveFatigue, FEElasticMaterial)
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
FEReactiveFatigue::FEReactiveFatigue(FEModel* pfem) : FEElasticMaterial(pfem)
{
	m_pBase = 0;
	m_pIdmg = 0;
	m_pFdmg = 0;
	m_pIcrt = 0;
	m_pFcrt = 0;
    
    m_k0 = 0;
    m_beta = 0;
}

//-----------------------------------------------------------------------------
// returns a pointer to a new material point object
FEMaterialPointData* FEReactiveFatigue::CreateMaterialPointData()
{
	return new FEReactiveFatigueMaterialPoint(m_pBase->CreateMaterialPointData());
}

//-----------------------------------------------------------------------------
//! Initialization.
bool FEReactiveFatigue::Init()
{
    FEUncoupledMaterial* m_pMat = dynamic_cast<FEUncoupledMaterial*>((FEElasticMaterial*)m_pBase);
	if (m_pMat != nullptr)
	{
		feLogError("Elastic material should not be of type uncoupled");
		return false;
	}
    
    return FEElasticMaterial::Init();
}

//-----------------------------------------------------------------------------
//! calculate stress at material point
mat3ds FEReactiveFatigue::Stress(FEMaterialPoint& pt)
{
    // evaluate the damage
    double d = Damage(pt);
    
    // evaluate the stress
    mat3ds s = m_pBase->Stress(pt);
    
    // return damaged stress
    return s*(1-d);
}

//-----------------------------------------------------------------------------
//! calculate tangent stiffness at material point
tens4ds FEReactiveFatigue::Tangent(FEMaterialPoint& pt)
{
    // evaluate the damage
    double d = Damage(pt);
    
    // evaluate the tangent
    tens4ds c = m_pBase->Tangent(pt);
    
    // return damaged tangent
    return c*(1-d);
}

//-----------------------------------------------------------------------------
//! calculate strain energy density at material point
double FEReactiveFatigue::StrainEnergyDensity(FEMaterialPoint& pt)
{
    // evaluate the damage
    double d = Damage(pt);
    
    // evaluate the strain energy density
    double sed = m_pBase->StrainEnergyDensity(pt);
    
    // return damaged sed
    return sed*(1-d);
}

//-----------------------------------------------------------------------------
//! calculate damage at material point
double FEReactiveFatigue::Damage(FEMaterialPoint& pt)
{
    // get the reactive fatigue material point data
    FEReactiveFatigueMaterialPoint& pd = *pt.ExtractData<FEReactiveFatigueMaterialPoint>();
    
    return pd.m_D;
}

//-----------------------------------------------------------------------------
// update fatigue material point at each iteration
void FEReactiveFatigue::UpdateSpecializedMaterialPoints(FEMaterialPoint& pt, const FETimeInfo& tp)
{
    double dt = tp.timeIncrement;
    double k0 = m_k0(pt);
    double beta = m_beta(pt);
    
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
    
    // evaluate increment in broken bond mass fraction from breaking of intact bonds
    double dwb_ib = (pd.m_Fip < 1) ? (pd.m_Fit-pd.m_Fip)/(1-pd.m_Fip)*pd.m_wip : pd.m_wip;
    // kinetics of intact bonds
    pd.m_wbt = pd.m_wbp + dwb_ib;
    // wit must be properly initialized before evaluating the fatigue reaction
    pd.m_wit = pd.m_wip - dwb_ib;

    // solve for fatigue bond mass fractions at current iteration
    double kt = k0*pow(fabs(pd.m_aXit)*pd.m_wbt,beta);
    double kp = k0*pow(fabs(pd.m_aXip)*pd.m_wbp,beta);
    // evaluate increment in fatigue bond mass fraction from fatigue of intact bonds
    double dwf_if = dt/2*(kt*pd.m_wit+kp*pd.m_wip);

    // kinetics of intact bonds
    pd.m_wit -= dwf_if;
    pd.m_wft = pd.m_wfp + dwf_if;
    
    // add or update new generation
    if ((pd.m_fb.size() == 0) || pd.m_fb.back().m_time < tp.currentTime) {
        // add generation of fatigued bonds
        FatigueBond fb;
        fb.m_Ffp = 0;
        fb.m_Fft = m_pFdmg->cdf(pt,Xftrl);
        fb.m_wfp = 0;
        fb.m_wft = dwf_if;
        fb.m_Xftrl = Xftrl;
        fb.m_time = tp.currentTime;
        pd.m_fb.push_back(fb);
    }
    else {
        for (int ig=0; ig < pd.m_fb.size(); ++ig) pd.m_fb[ig].m_wft = pd.m_fb[ig].m_wfp;
        pd.m_fb.back().m_Fft = m_pFdmg->cdf(pt,Xftrl);
        pd.m_fb.back().m_wft += dwf_if;
        pd.m_fb.back().m_Xftrl = Xftrl;
    }
    // damage kinetics of fatigued bonds
    for (int ig=0; ig < pd.m_fb.size(); ++ig) {
        double dwb_fb = (pd.m_fb[ig].m_Ffp < 1) ? (pd.m_fb[ig].m_Fft-pd.m_fb[ig].m_Ffp)/(1-pd.m_fb[ig].m_Ffp)*pd.m_fb[ig].m_wfp
        : pd.m_fb[ig].m_wft;
        pd.m_fb[ig].m_wft -= dwb_fb;
        pd.m_wbt += dwb_fb;
        pd.m_wft -= dwb_fb;
    }
    // evaluate fatigue bond fraction
//    pd.m_wft = 0;
//    for (int ig=0; ig < pd.m_fb.size(); ++ig) pd.m_wft += pd.m_fb[ig].m_wft;
    // fix round-off errors
    double wbd = 0;
    if (pd.m_wbt < pd.m_wbp) {
        wbd = pd.m_wbp - pd.m_wbt;
        pd.m_wbt = pd.m_wbp;
    }
    else if (pd.m_wbt > 1) {
        pd.m_wbt = 1;
    }
    double wid = 0;
    if (pd.m_wit > pd.m_wip) {
        wid = pd.m_wit - pd.m_wip;
        pd.m_wit = pd.m_wip;
    }
    else if (pd.m_wit < 0) {
        pd.m_wit = 0;
    }
    pd.m_wft = 1 - pd.m_wit - pd.m_wbt - wid + wbd;
    pd.m_D = pd.m_wbt;
}
