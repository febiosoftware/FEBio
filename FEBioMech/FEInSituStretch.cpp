#include "stdafx.h"
#include "FEInSituStretch.h"
#include "FEElasticMixture.h"
#include "FEUncoupledElasticMixture.h"
#include "FECore/FEModel.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_PARAMETER_LIST(FEInSituStretch, FENLConstraint)
	ADD_PARAMETER(m_ltol, FE_PARAM_DOUBLE, "augtol");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEInSituStretch::FEInSituStretch(FEModel* pfem) : FENLConstraint(pfem)
{
}

//-----------------------------------------------------------------------------
bool FEInSituStretch::Init()
{
	return true;
}

//-----------------------------------------------------------------------------
bool FEInSituStretch::Augment(int naug)
{
	if (IsActive() == false) return true;

	FEMesh& m = GetFEModel()->GetMesh();

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
			FEPreStrainTransIsoMR* pm = dynamic_cast<FEPreStrainTransIsoMR*>(pmat);
			if (pm) bconv |= Augment(psd, pm, 0);
			else
			{
				for (int j=0; j<pmat->Properties(); ++j)
				{
					FEPreStrainTransIsoMR* pmj = dynamic_cast<FEPreStrainTransIsoMR*>(pmat->GetProperty(j));
					if (pmj) bconv |= Augment(psd, pmj, j);
				}
			}
		}
	}

	return bconv;
}

//-----------------------------------------------------------------------------
//! Check an augmentation for a specific domain/material pair
bool FEInSituStretch::Augment(FESolidDomain* psd, FEPreStrainTransIsoMR* pmat, int n)
{
	int NE = psd->Elements();

	// check convergence
	bool bconv = true;
	for (int i=0; i<NE; ++i)
	{
		FESolidElement& el = psd->Element(i);
		int nint = el.GaussPoints();
		for (int i=0; i<nint; ++i)
		{
			FEMaterialPoint& mp = *(el.GetMaterialPoint(i)->GetPointData(n));
			FEFiberPreStretchMaterialPoint& psp = *mp.ExtractData<FEFiberPreStretchMaterialPoint>();

			FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
			mat3d& F = pt.m_F;

			vec3d a0(pt.m_Q[0][0], pt.m_Q[1][0], pt.m_Q[2][0]);
			vec3d a = F*a0;
			double lRtor = a.norm();

			// get the target stretch
			double ltrg = pmat->FiberStretch(mp);

			double err = ltrg / lRtor - psp.m_lam;
			double lnew = psp.m_lam + err;

			double dl = (lnew - psp.m_lamp)/psp.m_lamp;
			if (fabs(dl) >= m_ltol) bconv = false;
		}
	}

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
				FEFiberPreStretchMaterialPoint& psp = *mp.ExtractData<FEFiberPreStretchMaterialPoint>();

				FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
				mat3d& F = pt.m_F;

				vec3d a0(pt.m_Q[0][0], pt.m_Q[1][0], pt.m_Q[2][0]);
				vec3d a = F*a0;
				double lRtor = a.norm();

				// get the target stretch
				double ltrg = pmat->FiberStretch(mp);

				double err = ltrg / lRtor - psp.m_lam;
				double lnew = psp.m_lam + err;

				psp.m_lamp = psp.m_lam;
				psp.m_lam = lnew;
			}
		}
	}

	return bconv;
}
