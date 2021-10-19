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
#include <limits>
#include "FEFiberIntegrationGaussKronrod.h"
#include "gausskronrod.h"
#include <FECore/log.h>

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

class FEFiberIntegrationGaussKronrod::Iterator : public FEFiberIntegrationSchemeIterator
{
public:
	Iterator(FEMaterialPoint* mp, FEFiberIntegrationGaussKronrod::GKRULE& rule)
	{
		m_ncase = -1;
		m_nth = rule.m_nth;
		m_nph = rule.m_nph;
		m_gp = rule.m_gp;
		m_gw = rule.m_gw;

		const double eps = std::numeric_limits<double>::epsilon();
		if (mp)
		{
			FEElasticMaterialPoint& pt = *mp->ExtractData<FEElasticMaterialPoint>();

			// right Cauchy-Green tensor and its eigenvalues & eigenvectors
			// TODO: for uncoupled formulations we need to use the deviatoric strain-energy
			mat3ds C = (pt.m_buncoupled? pt.DevRightCauchyGreen() : pt.RightCauchyGreen());
			mat3ds E = (C - mat3dd(1)) * 0.5;
			E.eigen2(lE, vE);//lE[2]>lE[1]>lE[0]

			// check if there is no tension
			if (lE[2] <= 0) return;

			// check the other case
			if      (lE[0] >=   0) m_ncase = 0;
			else if (lE[1] <= eps) m_ncase = 1;
			else m_ncase = 2;
		}
		else
		{
			m_ncase = 0;
			vE[0] = vec3d(1,0,0);
			vE[1] = vec3d(0,1,0);
			vE[2] = vec3d(0,0,1);
		}

		// initialize theta increment
		double pi = 4 * atan(1.0);
		dth = 2 * pi / m_nth;

		// tension along all three eigenvectors (all directions)
		if (m_ncase == 0)
		{
			ksia = 0;
			ksib = 1;
			dksi = (ksib - ksia) / 2;
			sksi = (ksib + ksia) / 2;
		}
		// tension along one eigenvector and compression along other two// TC case
		else if (m_ncase == 1)
		{
			// nothing to do
			// ksia, ksib are calculated in FiberVector
		}
		// tension along two eigenvectors and compression along third
		else
		{
			// swap first and last eigenvalues/eigenvectors to maintain consistency in formulas
			double ltmp = lE[2]; vec3d vtmp = vE[2];
			lE[2] = lE[0]; vE[2] = vE[0];
			lE[0] = ltmp; vE[0] = vtmp;
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
		// check if the iterator is valid
		if (m_ncase == -1) return false;

		// update loop counters
		j++;
		if (j>=m_nph)
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
				double nu = -lE[0] * SQR(cth) - lE[1] * SQR(sth);
				if (nu<0) nu = 0;
				double de = lE[2] - lE[0] * SQR(cth) - lE[1] * SQR(sth);
				if (de <= 0) ksia = 1;
				else ksia = sqrt(nu) / sqrt(de);
				ksib = 1;
				dksi = (ksib - ksia) / 2;
				sksi = (ksib + ksia) / 2;
			}
			else if (m_ncase == 2)
			{
				ksia = 0;
				double nu = lE[0] * SQR(cth) + lE[1] * SQR(sth);
				if (nu<0) nu = 0;
				double de = lE[0] * SQR(cth) + lE[1] * SQR(sth) - lE[2];
				if (de <= 0) ksib = 0;
				else ksib = sqrt(nu) / sqrt(de);
				dksi = (ksib - ksia) / 2;
				sksi = (ksib + ksia) / 2;
			}
		}

		double ksi = sksi + dksi*m_gp[j];
		double sph = sqrt(1.0 - ksi*ksi); // = sin(acos(ksi));
		m_fiber = vE[0] * (cth*sph) + vE[1] * (sth*sph) + vE[2] * ksi;

		// we multiply by two to add contribution from other half-sphere
		m_weight = (m_gw[j] * dth*dksi)*2.0;

		i_old = i;
		return true;
	}

public:
	int	m_ncase;
	int	m_nth, m_nph;
	double lE[3];
	vec3d vE[3];

	int i, j, i_old;	// loop iterators
	double dth;
	double cth, sth;
	double ksia, ksib, dksi, sksi;
	const double*	m_gp;
	const double*	m_gw;
};

//-----------------------------------------------------------------------------
// FEFiberIntegrationGaussKronrod
//-----------------------------------------------------------------------------

// define the material parameters
BEGIN_FECORE_CLASS(FEFiberIntegrationGaussKronrod, FEFiberIntegrationScheme)
	ADD_PARAMETER(m_rule.m_nph, "nph");
	ADD_PARAMETER(m_rule.m_nth, "nth");
END_FECORE_CLASS();

FEFiberIntegrationGaussKronrod::FEFiberIntegrationGaussKronrod(FEModel* pfem) : FEFiberIntegrationScheme(pfem) 
{ 
	m_rule.m_nph = 7;
	m_rule.m_nth = 31;
}

//-----------------------------------------------------------------------------
FEFiberIntegrationGaussKronrod::~FEFiberIntegrationGaussKronrod()
{

}

//-----------------------------------------------------------------------------
void FEFiberIntegrationGaussKronrod::Serialize(DumpStream& ar)
{
	FEFiberIntegrationScheme::Serialize(ar);
	if ((ar.IsShallow() == false) && (ar.IsSaving() == false))
	{
		InitRule();
	}
}

//-----------------------------------------------------------------------------
bool FEFiberIntegrationGaussKronrod::Init()
{
	if (m_rule.m_nth < 1) {
		feLogError("nth must be strictly greater than zero."); return false;
	}

	// initialize the rule
	if (InitRule() == false) {
		feLogError("nph must be 7, 11, 15, 19, 23, or 27."); return false;
	}    
    
    // also initialize the parent class
    return FEFiberIntegrationScheme::Init();
}

//-----------------------------------------------------------------------------
bool FEFiberIntegrationGaussKronrod::InitRule()
{
	switch (m_rule.m_nph) {
	case 7:
		m_rule.m_gp = gp7;
		m_rule.m_gw = gw7;
		break;
	case 11:
		m_rule.m_gp = gp11;
		m_rule.m_gw = gw11;
		break;
	case 15:
		m_rule.m_gp = gp15;
		m_rule.m_gw = gw15;
		break;
	case 19:
		m_rule.m_gp = gp19;
		m_rule.m_gw = gw19;
		break;
	case 23:
		m_rule.m_gp = gp23;
		m_rule.m_gw = gw23;
		break;
	case 27:
		m_rule.m_gp = gp27;
		m_rule.m_gw = gw27;
		break;
	default:
		return false;
		break;
	}

	return true;
}

//-----------------------------------------------------------------------------
FEFiberIntegrationSchemeIterator* FEFiberIntegrationGaussKronrod::GetIterator(FEMaterialPoint* mp)
{
	// create a new iterator
	return new Iterator(mp, m_rule);
}
