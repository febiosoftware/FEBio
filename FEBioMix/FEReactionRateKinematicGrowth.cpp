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
#include "FEReactionRateKinematicGrowth.h"
#include "FEBiphasic.h"
#include <FEBioMech/FERemodelingElasticMaterial.h>
#include "FEReaction.h"
#include "FESoluteInterface.h"
#include <FECore\log.h>
#include <FEBioMech/FEGrowthTensor.h>
#include <FEBioMech/FEKinematicGrowth.h>
#include <FECore/FEModel.h>
#include <FECore/FEConstValueVec3.h>

// Material parameters for the FEMultiphasic material
BEGIN_FECORE_CLASS(FEReactionRateKinematicGrowth, FEReactionRate)
	ADD_PARAMETER(m_k, "k");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEReactionRateKinematicGrowth::FEReactionRateKinematicGrowth(FEModel* pfem) : FEReactionRate(pfem)
{

}

bool FEReactionRateKinematicGrowth::Init()
{
	/*if (m_fiber == nullptr) {
		FEConstValueVec3* val = fecore_new<FEConstValueVec3>("vector", nullptr);
		val->value() = vec3d(1, 0, 0);
		m_fiber = val;
	}*/

	return FEReactionRate::Init();
}

//-----------------------------------------------------------------------------
//! reaction rate at material point
double FEReactionRateKinematicGrowth::ReactionRate(FEMaterialPoint& pt)
{
	//FEDomain& dom = dynamic_cast<FEDomain&>(*pt.m_elem->GetMeshPartition());

	//FEModel* fem = this->m_pReact->GetFEModel();
	//FEMaterial* mat = dom.GetMaterial();
	//FEKinematicGrowth* kg = mat->ExtractProperty<FEKinematicGrowth>();
	//// Check if we were successful.
	//if (kg == nullptr) return false;
	//FEGrowthTensor* gmat = kg->GetGrowthMaterial();
	////FEVolumeGrowth* vg = dynamic_cast<FEVolumeGrowth*>(pmf);
	//if (gmat == nullptr) return false;
	////if (vg == nullptr) return false;
	//// For each element get the growth tensor and then solve the average value.
	//mat3d Q = gmat->GetLocalCS(pt);
	//vec3d fiber = gmat->m_fiber->unitVector(pt);
	//vec3d a0 = Q * fiber;
	//mat3d Fg = gmat->GrowthTensor(pt,a0);

	//mat3d Fgi = gmat->GrowthTensorInverse(pt, a0);
	//double Jgi = Fgi.det();

	// Get the deformation gradient and evaluate elastic deformation
	FEElasticMaterialPoint& ep = *pt.ExtractData<FEElasticMaterialPoint>();

	// Get the deformation gradient and evaluate elastic deformation
	FEKinematicMaterialPoint& kp = *pt.ExtractData<FEKinematicMaterialPoint>();

	FEBiphasicInterface* pbm = dynamic_cast<FEBiphasicInterface*>(GetAncestor());
	double phir = pbm->SolidReferentialVolumeFraction(pt);
	//SL: Workaround to get the actual concentration since every other function assumes c depends on Jtot rather than Je
	double zhat = m_k(pt) * (ep.m_J) / (kp.m_Je - phir);
	return zhat;
}

//-----------------------------------------------------------------------------
//! tangent of reaction rate with strain at material point
mat3ds FEReactionRateKinematicGrowth::Tangent_ReactionRate_Strain(FEMaterialPoint& pt)
{
	FEBiphasicInterface* pbm = dynamic_cast<FEBiphasicInterface*>(GetAncestor());
	double phir = pbm->SolidReferentialVolumeFraction(pt);
	FEElasticMaterialPoint& et = *pt.ExtractData<FEElasticMaterialPoint>();
	double J = et.m_J;
	vec3d pos = pt.m_r0;

	mat3ds I = mat3ds(1.0);
	double zhat = ReactionRate(pt);
	mat3ds dzde = -1.0 * zhat * I;
	//mat3ds dzde = -1.0 * (zhat / (J - phir)) * I;
	return dzde;
}

//-----------------------------------------------------------------------------
//! tangent of reaction rate with effective fluid pressure at material point
double FEReactionRateKinematicGrowth::Tangent_ReactionRate_Pressure(FEMaterialPoint& pt)
{
	return 0;
}

