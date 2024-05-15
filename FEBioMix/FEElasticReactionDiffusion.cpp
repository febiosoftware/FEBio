/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2024 University of Utah, The Trustees of Columbia University in
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
#include "FEElasticReactionDiffusion.h"
#include <FECore/FEModel.h>
#include <FECore/FECoreKernel.h>
#include <FECore/log.h>
#include <FECore/tens4d.h>
#include <FECore/tools.h>
#include <complex>
#include <FEBioMech/FEKinematicGrowth.h>
using namespace std;

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

// Material parameters for the FEElasticReactionDiffusion material
BEGIN_FECORE_CLASS(FEElasticReactionDiffusion, FEMaterial)
	ADD_PARAMETER(m_phi0   , FE_RANGE_CLOSED     (0.0, 1.0), "phi0"         );

	// define the material properties
	ADD_PROPERTY(m_pSolid , "solid"              , FEProperty::Required | FEProperty::TopLevel);
	ADD_PROPERTY(m_pSolute, "solute"             , FEProperty::Optional);
	ADD_PROPERTY(m_pReact , "reaction"           , FEProperty::Optional);

	ADD_PROPERTY(m_Q, "mat_axis", FEProperty::Optional);

END_FECORE_CLASS();

//=============================================================================
//   FEElasticReactionDiffusion
//=============================================================================

//-----------------------------------------------------------------------------
//! FEElasticReactionDiffusion constructor
FEElasticReactionDiffusion::FEElasticReactionDiffusion(FEModel* pfem) : FEMaterial(pfem)
{	
	m_Rgas = 0; m_Tabs = 0;

	m_pSolid = 0;
}

//-----------------------------------------------------------------------------
void FEElasticReactionDiffusion::AddChemicalReaction(FEChemicalReactionERD* pcr)
{
	m_pReact.push_back(pcr);
}

//-----------------------------------------------------------------------------
bool FEElasticReactionDiffusion::Init()
{
	// set the solute IDs first, since they are referenced in FESolute::Init()
	for (int i = 0; i<Solutes(); ++i) {
		m_pSolute[i]->SetSoluteLocalID(i);
	}

    if (m_pSolid->Init() == false) return false;
    for (int i=0; i<Solutes(); ++i)
        if (m_pSolute[i]->Init() == false) return false;
    for (int i=0; i<Reactions(); ++i)
        if (m_pReact[i]->Init() == false) return false;
    
	// call the base class.
	// This also initializes all properties
	if (FEMaterial::Init() == false) return false;

	m_Rgas = GetFEModel()->GetGlobalConstant("R");
	m_Tabs = GetFEModel()->GetGlobalConstant("T");

	if (m_Rgas <= 0) { feLogError("A positive universal gas constant R must be defined in Globals section"); return false; }
	if (m_Tabs <= 0) { feLogError("A positive absolute temperature T must be defined in Globals section");	 return false; }

	return true;
}

//-----------------------------------------------------------------------------
// update specialized material points
void FEElasticReactionDiffusion::UpdateSpecializedMaterialPoints(FEMaterialPoint& mp, const FETimeInfo& tp)
{
    m_pSolid->UpdateSpecializedMaterialPoints(mp, tp);
    for (int i=0; i<Solutes(); ++i)
        m_pSolute[i]->UpdateSpecializedMaterialPoints(mp, tp);
    for (int i=0; i<Reactions(); ++i)
        m_pReact[i]->UpdateSpecializedMaterialPoints(mp, tp);
}

//-----------------------------------------------------------------------------
void FEElasticReactionDiffusion::Serialize(DumpStream& ar)
{
	FEMaterial::Serialize(ar);
	if (ar.IsShallow()) return;

	ar & m_Rgas & m_Tabs;
}

//-----------------------------------------------------------------------------
//! Porosity in current configuration
double FEElasticReactionDiffusion::Porosity(FEMaterialPoint& pt)
{	
	// solid referential volume fraction
	double phisr = m_phi0(pt);
	
	double phiw = 1.0 - phisr;
	// check for pore collapse
	// TODO: throw an error if pores collapse
	phiw = (phiw > 0) ? phiw : 0;
	
	return phiw;
}

//-----------------------------------------------------------------------------
//! actual concentration
double FEElasticReactionDiffusion::Concentration(FEMaterialPoint& pt, const int sol)
{
	FESolutesMaterialPoint& spt = *pt.ExtractData<FESolutesMaterialPoint>();
	
	// effective and actual concentration
	double c = spt.m_c[sol];
	
	return c;
}

//-----------------------------------------------------------------------------
//! The stress of a triphasic material is the sum of the fluid pressure
//! and the elastic stress. Note that this function is declared in the base class
//! so you do not have to reimplement it in a derived class, unless additional
//! pressure terms are required.

mat3ds FEElasticReactionDiffusion::Stress(FEMaterialPoint& mp)
{
	// calculate solid material stress
	return m_pSolid->Stress(mp);
}

//-----------------------------------------------------------------------------
//! The tangent is the elastic tangent. Note
//! that this function is declared in the base class, so you don't have to 
//! reimplement it unless additional tangent components are required.

tens4ds FEElasticReactionDiffusion::Tangent(FEMaterialPoint& mp)
{	
	// call solid tangent routine
	return m_pSolid->Tangent(mp);
}

//-----------------------------------------------------------------------------
//! Calculate solute molar flux

vec3d FEElasticReactionDiffusion::SoluteFlux(FEMaterialPoint& pt, const int sol)
{
	FESolutesMaterialPoint& spt = *pt.ExtractData<FESolutesMaterialPoint>();
	
	// fluid volume fraction (porosity) in current configuration
	double phiw = Porosity(pt);
	
	// concentration
	double c = spt.m_c[sol];
	
	// concentration gradient
	vec3d gradc = spt.m_gradc[sol];
	
	// solute diffusivity in mixture
	mat3ds D = m_pSolute[sol]->m_pDiff->Diffusivity(pt);
	
	// solute flux j
    vec3d j = -D * gradc * phiw;
	
	return j;
}