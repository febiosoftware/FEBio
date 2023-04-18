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
#include "FEPreStrainConstraint.h"
#include "FEElasticMixture.h"
#include "FEUncoupledElasticMixture.h"
#include "FEPreStrainElastic.h"
#include <FECore/FECoreKernel.h>
#include <FECore/FEMesh.h>
#include <FECore/log.h>

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEPreStrainConstraint, FENLConstraint)
	ADD_PARAMETER(m_laugon , "update"   );
	ADD_PARAMETER(m_tol    , "tolerance");
	ADD_PARAMETER(m_naugmin, "min_iters");
	ADD_PARAMETER(m_naugmax, "max_iters");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEPreStrainConstraint::FEPreStrainConstraint(FEModel* pfem) : FENLConstraint(pfem)
{
	m_laugon = true; // This feature requires augmentations so turn it on
	m_naugmin = 0;
	m_naugmax = -1;
	m_tol = 0.0;
}

//-----------------------------------------------------------------------------
bool FEPreStrainConstraint::Init()
{
	return true;
}

//-----------------------------------------------------------------------------
bool FEPreStrainConstraint::Augment(int naug, const FETimeInfo& tp)
{
	if (IsActive() == false) return true;
	if (m_laugon == false) return true;

	FEMesh& m = GetMesh();
	int ND = m.Domains();

	// do pre-strain augmentations
	bool bconv = true;
	for (int i=0; i<ND; ++i)
	{
		FEDomain& dom = m.Domain(i);
		if (dom.Class() == FE_DOMAIN_SOLID)
		{
			FESolidDomain* psd = static_cast<FESolidDomain*>(&dom);
			FEMaterial* pmat = psd->GetMaterial();
			FEPrestrainMaterial* pm = dynamic_cast<FEPrestrainMaterial*>(pmat);
			if (pm) bconv &= Augment(psd, 0, naug);
			else
			{
				for (int j=0; j<pmat->Properties(); ++j)
				{
					FEPrestrainMaterial* pmj = dynamic_cast<FEPrestrainMaterial*>(pmat->GetProperty(j));
					if (pmj) bconv &= Augment(psd, j, naug);
				}
			}
		}
	}

	return bconv;
}

//-----------------------------------------------------------------------------
//! Check an augmentation for a specific domain/material pair
bool FEPreStrainConstraint::Augment(FESolidDomain* psd, int n, int naug)
{
	// make sure this is a prestrain material
	FEPrestrainMaterial* pmat = dynamic_cast<FEPrestrainMaterial*>(psd->GetMaterial());
	if (pmat == nullptr) return true;

	// check convergence
	double	max_err = 0;
	bool bconv = true;
	int NE = psd->Elements();
	for (int i=0; i<NE; ++i)
	{
		FESolidElement& el = psd->Element(i);
		int nint = el.GaussPoints();
		for (int i=0; i<nint; ++i)
		{
			FEMaterialPoint& mp = *(el.GetMaterialPoint(i)->GetPointData(n));
			FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();
			FEPrestrainMaterialPoint& pt = *mp.ExtractData<FEPrestrainMaterialPoint>();

			const mat3d& Fc = pt.PrestrainCorrection();
			const mat3d& F = ep.m_F;

			mat3d Fc_next = UpdateFc(F, Fc, mp, pmat);

			mat3d U = Fc_next - Fc;
			double normU = U.norm();
			if (normU >= max_err) max_err = normU;
		}
	}
	if (max_err >= m_tol) bconv = false;

	feLog("max norm = %lg (%lg)\n", max_err, m_tol);

	// ensure we have done the required min or max augmentations
	if (naug <  m_naugmin) bconv = false;
	if ((m_naugmax >= 0)&&(naug >= m_naugmax)) bconv = true;

	// only augment when we did not converge
	if (bconv == false)
	{
		for (int i=0; i<NE; ++i)
		{
			FESolidElement& el = psd->Element(i);
			int nint = el.GaussPoints();
			for (int i=0; i<nint; ++i)
			{
				FEMaterialPoint& mp = *(el.GetMaterialPoint(i)->GetPointData(n));
				FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();
				FEPrestrainMaterialPoint& pt = *mp.ExtractData<FEPrestrainMaterialPoint>();

				const mat3d& Fc = pt.PrestrainCorrection();
				const mat3d& F = ep.m_F;

				mat3d Fc_next = UpdateFc(F, Fc, mp, pmat);
				pt.setPrestrainCorrection(Fc_next);
			}
		}
	}

	return bconv;
}

//-----------------------------------------------------------------------------
mat3d FEGPAConstraint::UpdateFc(const mat3d& F, const mat3d& Fc_prev, FEMaterialPoint& mp, FEPrestrainMaterial* pmat)
{
	return F * Fc_prev;
}

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEInSituStretchConstraint, FEPreStrainConstraint)
	ADD_PARAMETER(m_max_stretch, "max_stretch");
	ADD_PARAMETER(m_biso       , "isochoric"  );
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEInSituStretchConstraint::FEInSituStretchConstraint(FEModel* pfem) : FEPreStrainConstraint(pfem)
{
	m_max_stretch = 0.0;
	m_biso = true;
}

//-----------------------------------------------------------------------------
mat3d FEInSituStretchConstraint::UpdateFc(const mat3d& F, const mat3d& Fc_prev, FEMaterialPoint& mp, FEPrestrainMaterial* pmat)
{
	FEPrestrainMaterialPoint& psp = *mp.ExtractData<FEPrestrainMaterialPoint>();
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	FEElasticMaterial* elasticMat = pmat->GetElasticMaterial();
	if (elasticMat == nullptr) return mat3dd(1.0);

	FEVec3dValuator* fiber = dynamic_cast<FEVec3dValuator*>(elasticMat->GetProperty("fiber"));
	if (fiber == nullptr) return mat3dd(1.0);

	// calculate the fiber stretch
	mat3d Q = elasticMat->GetLocalCS(mp);
	vec3d a0 = fiber->unitVector(mp);
	vec3d ar = Q * a0;
	vec3d a = F*ar;
	double l = a.norm();

	if ((m_max_stretch == 0.0) ||(l <= m_max_stretch))
	{
		double li = (m_biso ? sqrt(l) : 1.0);

		// setup the new correction
		mat3d U(1.0/l, 0.0, 0.0, 0.0, li, 0.0, 0.0, 0.0, li);
		mat3d Fc = Q*U*Q.transpose();
		return Fc;
	}
	else
	{
		return psp.PrestrainCorrection();
	}
}
