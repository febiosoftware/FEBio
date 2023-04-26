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
#include "FEReactionRateStressSensitive.h"
#include "FEBiphasic.h"
#include "FEBioMech/FERemodelingElasticMaterial.h"
#include <FECore\log.h>
#include "FEBioMech/FENeoHookean.h"
#include <FEBioMech/FEKinematicGrowth.h>
#include <FEBioMech/FEGrowthTensor.h>
#include <FECore\FESolidElement.h>
#include <FECore\FESPRProjection.h>
#include <FECore\FESolidDomain.h>

// Material parameters for the FEMultiphasic material
BEGIN_FECORE_CLASS(FEReactionRateStressSensitive, FEReactionRate)
	ADD_PARAMETER(m_k, "k");
	ADD_PARAMETER(m_a0, "a0");
	ADD_PARAMETER(m_a, "a");
	ADD_PARAMETER(m_b, "b");
	ADD_PARAMETER(stress0, "residual_stress")->setLongName("initial residual stress");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEReactionRateStressSensitive::FEReactionRateStressSensitive(FEModel* pfem) : FEReactionRate(pfem)
{ 

}

//-----------------------------------------------------------------------------
//! reaction rate at material point
double FEReactionRateStressSensitive::ReactionRate(FEMaterialPoint& pt)
{
	FEElasticMaterialPoint& ep = *pt.ExtractData<FEElasticMaterialPoint>();
	double m_S = m_a0 + m_a / (1.0 + (exp(-(ep.m_s.tr() - m_b) / stress0)));
	double zhat = m_k(pt) * m_S;
	return zhat;
}

//-----------------------------------------------------------------------------
//! tangent of reaction rate with strain at material point
//! SL! Todo: Figure out what to do with rhor for solutes.
mat3ds FEReactionRateStressSensitive::GetSPRStress(FEMaterialPoint& pt)
{
	const int LUT[6][2] = { { 0,0 },{ 1,1 },{ 2,2 },{ 0,1 },{ 1,2 },{ 0,2 } };
	
	int NN = pt.m_elem->Nodes();
	int NE = 1;

	FESolidElement& el = dynamic_cast<FESolidElement&>(*pt.m_elem);

	// build the element data array
	vector< vector<double> > ED[6];
	for (int n = 0; n < 6; ++n)
	{
		ED[n].resize(NE);

		int nint = el.GaussPoints();
		ED[n][0].assign(nint, 0.0);
	}

	// this array will store the results
	FESPRProjection map;
	map.SetInterpolationOrder(2);
	vector<double> val[6];

	int nint = el.GaussPoints();
	for (int i = 0; i < nint; ++i)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(i);
		FEElasticMaterialPoint& ep = dynamic_cast<FEElasticMaterialPoint&>(mp);
		mat3ds s = ep.m_s;

		// loop over stress components
		for (int n = 0; n < 6; ++n)
		{
			ED[n][0][i] = s(LUT[n][0], LUT[n][1]);
		}
	}
	
	// project to nodes
	// loop over stress components
	FEDomain& dom = dynamic_cast<FEDomain&>(*el.GetMeshPartition());
	FESolidDomain& sd = static_cast<FESolidDomain&>(dom);
	for (int n = 0; n < 6; ++n)
	{
		map.Project(sd, ED[n], val[n]);
	}
	mat3ds spr_s;
	// copy results to archive
	//for (int i = 0; i < NN; ++i)
	//{
	//	spr_s() = val[0][i];
	//	spr.push_back((float)val[0][i]);
	//	spr.push_back((float)val[1][i]);
	//	spr.push_back((float)val[2][i]);
	//	spr.push_back((float)val[3][i]);
	//	spr.push_back((float)val[4][i]);
	//	spr.push_back((float)val[5][i]);
	//}
	return mat3ds(0.0);
}

//-----------------------------------------------------------------------------
//! tangent of reaction rate with strain at material point
//! SL! Todo: Figure out what to do with rhor for solutes.
mat3ds FEReactionRateStressSensitive::Tangent_ReactionRate_Strain(FEMaterialPoint& pt)
{
	FEBiphasicInterface* pbm = dynamic_cast<FEBiphasicInterface*>(GetAncestor());
	double phir = pbm->SolidReferentialVolumeFraction(pt);
	FEElasticMaterialPoint& et = *pt.ExtractData<FEElasticMaterialPoint>();
	double J = et.m_J;
	vec3d pos = pt.m_r0;

	mat3ds I = mat3ds(1.0);
	double zhat = ReactionRate(pt);
	mat3ds dzde = -1.0 * zhat * I;
	return dzde;
}

//-----------------------------------------------------------------------------
//! tangent of reaction rate with effective fluid pressure at material point
double FEReactionRateStressSensitive::Tangent_ReactionRate_Pressure(FEMaterialPoint& pt)
{
	return 0;
}

