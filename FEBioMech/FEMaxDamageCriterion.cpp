/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2020 University of Utah, The Trustees of Columbia University in 
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
#include "FEMaxDamageCriterion.h"
#include "FEElasticMaterial.h"
#include <FECore/FEElement.h>
#include <FECore/FEModel.h>

BEGIN_FECORE_CLASS(FEMaxDamageCriterion, FEMeshAdaptorCriterion)
	ADD_PARAMETER(m_maxDamage, FE_RANGE_GREATER(0.0), "max_damage");
END_FECORE_CLASS();

FEMaxDamageCriterion::FEMaxDamageCriterion(FEModel* fem) : FEMeshAdaptorCriterion(fem)
{
	m_maxDamage = 0.0;

	// set sort on by default
	SetSort(true);
}

bool FEMaxDamageCriterion::Check(FEElement& el, double& elemVal)
{
	// return true if any of the integration point value exceeds the threshold
	bool bselect = false;
	int nint = el.GaussPoints();
	elemVal = 0;
	for (int n = 0; n < nint; ++n)
	{
		FEDamageMaterialPoint* dp = nullptr;
		FEMaterialPoint* mp = el.GetMaterialPoint(n);
		// for mixtures, we have to make sure we get the right component
		if (mp->Components() > 1)
		{
			for (int i = 0; i < mp->Components(); ++i)
			{
				FEMaterialPoint* mpi = mp->GetPointData(i);
				dp = mpi->ExtractData<FEDamageMaterialPoint>();
				if (dp) break;
			}
		}
		else dp = mp->ExtractData<FEDamageMaterialPoint>();

		// evaluate the damage at this point
		if (dp)
		{
			elemVal = dp->m_D;

			if (elemVal >= m_maxDamage)
			{
				bselect = true;
				break;
			}
		}
	}

	return bselect;
}
