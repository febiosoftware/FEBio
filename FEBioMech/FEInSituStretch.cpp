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
		FESolidDomain* psd = dynamic_cast<FESolidDomain*>(&dom);
		if (psd)
		{
			FEPreStrainTransIsoMR* pm = dynamic_cast<FEPreStrainTransIsoMR*>(psd->GetMaterial());
			if (pm) bconv = CheckAugment(psd, pm, 0);
			else
			{
				FEElasticMixture* pmix = dynamic_cast<FEElasticMixture*>(psd->GetMaterial());
				if (pmix)
				{
					for (int j=0; j<pmix->Materials(); ++j)
					{
						FEPreStrainTransIsoMR* pmj = dynamic_cast<FEPreStrainTransIsoMR*>(pmix->GetMaterial(j));
						if (pmj) bconv = CheckAugment(psd, pmj, j);
					}
				}
				else
				{
					FEUncoupledElasticMixture* pmix = dynamic_cast<FEUncoupledElasticMixture*>(psd->GetMaterial());
					if (pmix)
					{
						for (int j=0; j<pmix->Materials(); ++j)
						{
							FEPreStrainTransIsoMR* pmj = dynamic_cast<FEPreStrainTransIsoMR*>(pmix->GetMaterial(j));
							if (pmj) bconv = CheckAugment(psd, pmj, j);
						}
					}
				}
			}
		}
	}

	// if not converged, update Lagrange multipliers
	if (bconv == false)
	{
		for (int i=0; i<ND; ++i)
		{
			FEDomain& dom = m.Domain(i);
			FESolidDomain* psd = dynamic_cast<FESolidDomain*>(&dom);
			if (psd)
			{
				FEPreStrainTransIsoMR* pm = dynamic_cast<FEPreStrainTransIsoMR*>(psd->GetMaterial());
				if (pm) DoAugment(psd, pm, 0);
				else
				{
					FEElasticMixture* pmix = dynamic_cast<FEElasticMixture*>(psd->GetMaterial());
					if (pmix)
					{
						for (int j=0; j<pmix->Materials(); ++j)
						{
							FEPreStrainTransIsoMR* pmj = dynamic_cast<FEPreStrainTransIsoMR*>(pmix->GetMaterial(j));
							if (pmj) DoAugment(psd, pmj, j);
						}
					}
					else
					{
						FEUncoupledElasticMixture* pmix = dynamic_cast<FEUncoupledElasticMixture*>(psd->GetMaterial());
						if (pmix)
						{
							for (int j=0; j<pmix->Materials(); ++j)
							{
								FEPreStrainTransIsoMR* pmj = dynamic_cast<FEPreStrainTransIsoMR*>(pmix->GetMaterial(j));
								if (pmj) DoAugment(psd, pmj, j);
							}
						}
					}
				}
			}
		}
	}

	return bconv;
}

//-----------------------------------------------------------------------------
//! Check an augmentation for a specific domain/material pair
bool FEInSituStretch::CheckAugment(FESolidDomain* psd, FEPreStrainTransIsoMR* pmat, int n)
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

	return bconv;
}

//-----------------------------------------------------------------------------
void FEInSituStretch::DoAugment(FESolidDomain* psd, FEPreStrainTransIsoMR* pmat, int n)
{
	int NE = psd->Elements();
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
