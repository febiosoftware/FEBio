/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, Columbia University, and others.

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
