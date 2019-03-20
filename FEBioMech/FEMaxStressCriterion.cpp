#include "stdafx.h"
#include "FEMaxStressCriterion.h"
#include "FEElasticMaterial.h"
#include <FECore/FEElement.h>

BEGIN_FECORE_CLASS(FEMaxStressCriterion, FEMeshAdaptorCriterion)
	ADD_PARAMETER(m_maxStress, FE_RANGE_GREATER(0.0), "max_stress");
	ADD_PARAMETER(m_metric, "metric");
END_FECORE_CLASS();

FEMaxStressCriterion::FEMaxStressCriterion(FEModel* fem) : FEMeshAdaptorCriterion(fem)
{
	m_maxStress = 0.0;
	m_metric = 0;

	// set sort on by default
	SetSort(true);
}

bool FEMaxStressCriterion::Check(FEElement& el, double& elemVal)
{
	bool bselect = false;
	int nint = el.GaussPoints();
	elemVal = 0;
	for (int n = 0; n < nint; ++n)
	{
		FEMaterialPoint* mp = el.GetMaterialPoint(n);
		FEElasticMaterialPoint* ep = mp->ExtractData<FEElasticMaterialPoint>();
		if (ep)
		{
			mat3ds& s = ep->m_s;

			switch (m_metric)
			{
			case 0: elemVal = s.effective_norm(); break;
			case 1:
			{
				mat3ds devs = s.dev();
				double l[3], lmax;
				devs.exact_eigen(l);
				lmax = l[0];
				if (l[1] > lmax) lmax = l[1];
				if (l[2] > lmax) lmax = l[2];
				elemVal = lmax;
			}
			break;
			}

			if (elemVal >= m_maxStress)
			{
				bselect = true;
				break;
			}
		}
	}

	return bselect;
}
