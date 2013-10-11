#include "stdafx.h"
#include "FEInSituStretch.h"
#include "FEPreStrainTransIsoMR.h"
#include "FECore/FEModel.h"
#include "FECore/FESolidDomain.h"

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
void FEInSituStretch::Init()
{
}

//-----------------------------------------------------------------------------
bool FEInSituStretch::Augment(int naug)
{
	FEMesh& m = m_pfem->GetMesh();

	int ND = m.Domains();

	// do pre-strain augmentations
	bool bconv = true;

	for (int i=0; i<ND; ++i)
	{
		FEDomain& dom = m.Domain(i);
		FESolidDomain* psd = dynamic_cast<FESolidDomain*>(&dom);
		FEPreStrainTransIsoMR* pm = dynamic_cast<FEPreStrainTransIsoMR*>(dom.GetMaterial());
		if (psd && pm)
		{
			int NE = psd->Elements();

			// check convergence
			for (int i=0; i<NE; ++i)
			{
				FESolidElement& el = psd->Element(i);
				int nint = el.GaussPoints();
				for (int i=0; i<nint; ++i)
				{
					FEMaterialPoint& mp = *el.m_State[i];
					FEPreStrainMaterialPoint& psp = *mp.ExtractData<FEPreStrainMaterialPoint>();

					FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
					mat3d& F = pt.m_F;

					vec3d a0(pt.m_Q[0][0], pt.m_Q[1][0], pt.m_Q[2][0]);
					vec3d a = F*a0;
					double lRtor = a.norm();

					double err = pm->m_ltrg / lRtor - psp.m_lam;
					double lnew = psp.m_lam + err;

					double dl = (lnew - psp.m_lamp)/psp.m_lamp;
					if (fabs(dl) >= m_ltol) bconv = false;
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
			FEPreStrainTransIsoMR* pm = dynamic_cast<FEPreStrainTransIsoMR*>(dom.GetMaterial());
			if (psd && pm)
			{
				int NE = psd->Elements();

				for (int i=0; i<NE; ++i)
				{
					FESolidElement& el = psd->Element(i);
					int nint = el.GaussPoints();
					for (int i=0; i<nint; ++i)
					{
						FEMaterialPoint& mp = *el.m_State[i];
						FEPreStrainMaterialPoint& psp = *mp.ExtractData<FEPreStrainMaterialPoint>();

						FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
						mat3d& F = pt.m_F;

						vec3d a0(pt.m_Q[0][0], pt.m_Q[1][0], pt.m_Q[2][0]);
						vec3d a = F*a0;
						double lRtor = a.norm();

						double err = pm->m_ltrg / lRtor - psp.m_lam;
						double lnew = psp.m_lam + err;

						psp.m_lamp = psp.m_lam;
						psp.m_lam = lnew;
					}
				}
			}
		}
	}

	return bconv;
}
