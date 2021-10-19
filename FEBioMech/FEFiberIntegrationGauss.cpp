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
#include "FEFiberIntegrationGauss.h"
#include "gauss.h"
#include <FECore/log.h>

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

class FEFiberIntegrationGauss::Iterator : public FEFiberIntegrationSchemeIterator
{
public:
	Iterator(FEMaterialPoint* mp, FEFiberIntegrationGauss::GRULE& rule)
	{
		m_ncase = -1;
		m_nth = rule.m_nth;
		m_nph = rule.m_nph;
		m_gp = rule.m_gp;
		m_gw = rule.m_gw;

		const double eps = 1e-9;
		if (mp)
		{
			FEElasticMaterialPoint& pt = *mp->ExtractData<FEElasticMaterialPoint>();
			// right Cauchy-Green tensor and its eigenvalues & eigenvectors
			mat3ds C = (pt.m_buncoupled ? pt.DevRightCauchyGreen() : pt.RightCauchyGreen());
			C.eigen(lC, vC);

			// check if there is no tension
			if ((lC[0] <= 1 + eps) && (lC[1] <= 1 + eps) && (lC[2] <= 1 + eps)) {
				return;
			}

			// bubble sort eigenvalues & eigenvectors from smallest to largest
			double ltmp;
			vec3d vtmp;
			bool swp = true;
			while (swp) {
				swp = false;
				for (int i = 1; i<3; ++i) {
					int j = i - 1;
					if (lC[i] < lC[j]) {
						ltmp = lC[i]; vtmp = vC[i];
						lC[i] = lC[j]; vC[i] = vC[j];
						lC[j] = ltmp; vC[j] = vtmp;
						swp = true;
					}
				}
			}

			// check the other case
			if      (lC[0]  > 1 + eps) m_ncase = 0;
			else if (lC[1] <= 1 + eps) m_ncase = 1;
			else m_ncase = 2;
		}
		else
		{
			m_ncase = 0;
			vC[0] = vec3d(1,0,0);
			vC[1] = vec3d(0,1,0);
			vC[2] = vec3d(0,0,1);
		}

		// check remaining stretch states
		double pi = 4 * atan(1.0);
		dth = 2 * pi / m_nth;

		vec3d n0e, n0a;

		// tension along all three eigenvectors (all directions)
		if (m_ncase == 0) 
		{
			ksia = 0;
			ksib = 1;
			dksi = (ksib - ksia) / 2;
			sksi = (ksib + ksia) / 2;
		}
		// tension along one eigenvector and compression along other two
		else if (m_ncase == 1)
		{
			// nothing to do
			// ksia and ksib are update in FiberVector
		}
		// tension along two eigenvectors and compression along third
		else
		{
			// swap first and last eigenvalues/eigenvectors to maintain consistency in formulas
			double ltmp = lC[2]; vec3d vtmp = vC[2];
			lC[2] = lC[0]; vC[2] = vC[0];
			lC[0] = ltmp; vC[0] = vtmp;
		}

		i = 0; j = -1;
		i_old = -1;
		cth = 1.0; sth = 0.0;
		Next();
	}

	bool IsValid()
	{
		return (m_ncase != -1);
	}

	// move to the next integration point
	bool Next()
	{
		// make sure the iterator is valid
		if (m_ncase == -1) return false;

		// update loop counters
		j++;
		if (j >= m_nph)
		{
			j = 0;
			i++;
			if (i >= m_nth)
			{
				// all done
				m_ncase = -1;
				return false;
			}
		}

		// update vector
		if (i_old != i)
		{
			double theta = i*dth;
			cth = cos(theta);
			sth = sin(theta);

			if (m_ncase == 1)
			{
				ksia = sqrt(1 - lC[0] * SQR(cth) - lC[1] * SQR(sth)) / sqrt(lC[2] - lC[0] * SQR(cth) - lC[1] * SQR(sth));
				ksib = 1;
				dksi = (ksib - ksia) / 2;
				sksi = (ksib + ksia) / 2;
			}
			else if (m_ncase == 2)
			{
				ksia = 0;
				ksib = sqrt(lC[0] * SQR(cth) + lC[1] * SQR(sth) - 1) / sqrt(lC[0] * SQR(cth) + lC[1] * SQR(sth) - lC[2]);
				dksi = (ksib - ksia) / 2;
				sksi = (ksib + ksia) / 2;
			}
		}

		double ksi = sksi + dksi*m_gp[j];
		double sph = sqrt(1.0 - ksi*ksi); // = sin(acos(ksi));

		m_fiber = vC[0] * (cth*sph) + vC[1] * (sth*sph) + vC[2] * ksi;

		// we multiply by two to add contribution from other half-sphere
		m_weight = (m_gw[j] * dth*dksi)*2.0;

		i_old = i;

		return true;
	}

public:
	int m_ncase;
	double lC[3];
	vec3d vC[3];
	int	m_nth, m_nph;
	const double* m_gp;
	const double* m_gw;

	int i, j, i_old;
	double ksia, ksib, dksi, sksi;
	double dth;
	double cth, sth;
};

//-----------------------------------------------------------------------------
// FEFiberIntegrationGauss
//-----------------------------------------------------------------------------

// define the material parameters
BEGIN_FECORE_CLASS(FEFiberIntegrationGauss, FEFiberIntegrationScheme)
	ADD_PARAMETER(m_rule.m_nph, FE_RANGE_GREATER(0), "nph");
	ADD_PARAMETER(m_rule.m_nth, FE_RANGE_GREATER(0), "nth");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
void FEFiberIntegrationGauss::Serialize(DumpStream& ar)
{
	FEFiberIntegrationScheme::Serialize(ar);
	if ((ar.IsSaving() == false) && (ar.IsShallow() == false))
	{
		InitRule();
	}
}

//-----------------------------------------------------------------------------
bool FEFiberIntegrationGauss::InitRule()
{
	switch (m_rule.m_nph) {
        case 1:
			m_rule.m_gp = gp1;
			m_rule.m_gw = gw1;
            break;
        case 2:
			m_rule.m_gp = gp2;
			m_rule.m_gw = gw2;
            break;
        case 3:
			m_rule.m_gp = gp3;
			m_rule.m_gw = gw3;
            break;
        case 4:
			m_rule.m_gp = gp4;
			m_rule.m_gw = gw4;
            break;
        case 5:
			m_rule.m_gp = gp5;
			m_rule.m_gw = gw5;
            break;
        case 6:
			m_rule.m_gp = gp6;
			m_rule.m_gw = gw6;
            break;
        case 7:
			m_rule.m_gp = gp7;
			m_rule.m_gw = gw7;
            break;
        case 8:
			m_rule.m_gp = gp8;
			m_rule.m_gw = gw8;
            break;
        case 9:
			m_rule.m_gp = gp9;
			m_rule.m_gw = gw9;
            break;
        case 10:
			m_rule.m_gp = gp10;
			m_rule.m_gw = gw10;
            break;
        default:
            return false;
            break;
    }

	return true;
}

FEFiberIntegrationGauss::FEFiberIntegrationGauss(FEModel* pfem) : FEFiberIntegrationScheme(pfem)
{ 
	m_rule.m_nph = 5; 
	m_rule.m_nth = 2 * m_rule.m_nph;
}

FEFiberIntegrationGauss::~FEFiberIntegrationGauss()
{
}

bool FEFiberIntegrationGauss::Init()
{
	if (InitRule() == false) {
		feLogError("nint must not exceed 10."); return false;
	}

    // also initialize the parent class
    return FEFiberIntegrationScheme::Init();
}

//-----------------------------------------------------------------------------
FEFiberIntegrationSchemeIterator* FEFiberIntegrationGauss::GetIterator(FEMaterialPoint* mp)
{
	return new Iterator(mp, m_rule);
}
