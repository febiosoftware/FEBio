/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in 
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
#include "FEFiberEntropyChain.h"
#include <iostream>
#include "triangle_sphere.h"
#include <limits>
#include <FECore/log.h>

//-----------------------------------------------------------------------------
// FEFiberEntropyChain
//-----------------------------------------------------------------------------

// define the material parameters
BEGIN_FECORE_CLASS(FEFiberEntropyChain, FEFiberMaterial)
	ADD_PARAMETER(m_N   , FE_RANGE_GREATER_OR_EQUAL(0.0), "N");
	ADD_PARAMETER(m_ksi, FE_RANGE_GREATER_OR_EQUAL(0.0), "ksi")->setUnits(UNIT_PRESSURE);
    ADD_PARAMETER(m_mu , FE_RANGE_GREATER_OR_EQUAL(0.0), "mu" )->setUnits(UNIT_PRESSURE);
	ADD_PARAMETER(m_term, FE_RANGE_CLOSED(3,30), "n_term");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEFiberEntropyChain::FEFiberEntropyChain(FEModel* pfem) : FEFiberMaterial(pfem)
{
	m_mu = 0;
	m_term = 30;

	m_epsf = 1.0;
}

//-----------------------------------------------------------------------------
bool FEFiberEntropyChain::Validate()
{
    return FEFiberMaterial::Validate();
}

//-----------------------------------------------------------------------------
mat3ds FEFiberEntropyChain::FiberStress(FEMaterialPoint& mp, const vec3d& n0)
{
	double ksi = m_ksi(mp);
	double mu = m_mu(mp);
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    	
	// deformation gradient
	mat3d &F = pt.m_F;
	double J = pt.m_J;
	
	// loop over all integration points
	mat3ds C = pt.RightCauchyGreen();
	mat3ds s;
	
	const double f_a[] = { 3.0, 9.0 / 5.0, 297.0 / 175.0, 1539.0 / 875.0, 126117.0 / 67373.0, 43733439.0 / 21896875.0, \
		231321177.0 / 109484375.0, 20495009043.0 / 9306171875.0, 1073585186448381.0 / 476522530859375.0, 4387445039583.0 / 1944989921875.0, \
		1000263375846831627.0 / 453346207767578125.0, 280865021365240713.0 / 133337119931640625.0, 148014872740758343473.0 / 75350125192138671875.0, 137372931237386537808993.0 / 76480377070020751953125.0, \
		41722474198742657618857192737.0 / 25674386102028896409912109375.0, 12348948373636682700768301125723.0 / 8344175483159391333221435546875.0, 5001286000741585238340074032091091.0 / 3590449627018291035442047119140625.0, \
		185364329915163811141785118512534489.0 / 132846636199676768311355743408203125.0, 6292216384025878939310787532157558451.0 / 4160197291516193533960877227783203125.0, 299869254759556271677902570230858640837.0 / 170568088952163934892395966339111328125.0, \
		316689568216860631885141475537178451746044283.0 / 148810301691076651811854156590061187744140625.0, 670194310437429598283653289122392937145371137.0 / 258140319260030926612400067554187774658203125.0, 19697015384373759058671314622426656486031717197919.0 / 6332687088594936951180328352888571262359619140625.0, \
		178793788985653424246012689916144867915861856840849.0 / 49756827124674504616416865629838774204254150390625.0, 323844166067349737493036492206152479344269351967043143667.0 / 82039106319990447886052744467343800746748447418212890625.0, 200808116689754604893460969866238617668631975356485302537199.0 / 49409916306357883385918130190559334540655314922332763671875.0, \
		27506481209689719715275759452624078040221544551995885750037973.0 / 7164437864421893090958128877631103508395020663738250732421875.0, 16356939619211770477227805130221533318985185730316048281126247721.0 / 5122573073061653560035062147506239008502439774572849273681640625.0, 30976199222209837888906735596203249520053107250807132769871859115868101.0 / 14801698741072505317694781203956860833766632412399612367153167724609375.0, \
		519588001407316958447129785511020819131555326399179970047767492196701159.0 / 902903623205422824379381653441368510859764577156376354396343231201171875.0 };

	// fiber direction in global coordinate system
	//vec3d n0 = GetFiberVector(mp);
	
	// Calculate In = n0*C*n0
	double In = n0*(C*n0);
	
	// only take fibers in tension into consideration
	const double eps = m_epsf * std::numeric_limits<double>::epsilon();

	if ((In - 1.0) > eps)
	{
		// define the distribution
		double R, Ra = 0.0;
		// get the global spatial fiber direction in current configuration
		vec3d nt = F*n0;
		
		// calculate the outer product of nt
		mat3ds N = dyad(nt);
		
		// calculate strain energy derivative
		//double alpha = sqrt(In) / sqrt(m_N);
		double alpha = sqrt(In) / sqrt(m_N);
		double alpha0 = 1 / sqrt(m_N);
		double alphaI = 1 / sqrt(In) / (2.0 * sqrt(m_N));

		///////////////////////////////////////////////////////////////////////////////////////
		double beta = 0, beta0 = 0, alpha2 = alpha * alpha, alpha02 = alpha0 * alpha0; 
		int num_fa = sizeof(f_a) / sizeof(f_a[0]);
		if (m_term > num_fa) {
			m_term = num_fa;
		}
		beta = 1 + f_a[m_term - 1] / f_a[m_term - 2] * alpha2;
		beta0 = 1 + f_a[m_term - 1] / f_a[m_term - 2] * alpha02;
		for (int i = 2; i < m_term; i++) {
			beta = 1 + f_a[m_term - i] / f_a[m_term - i - 1] * alpha2 * beta;
			beta0 = 1 + f_a[m_term - i] / f_a[m_term - i - 1] * alpha02 * beta0;
			//cout << num_fa - i << endl;
		};
		beta = f_a[0] * alpha * beta;
		beta0 = f_a[0] * alpha0 * beta0;
		////////////////////////////////////////////////////////////////////////////////////////////
		
		double Wl = ksi * m_N * alphaI * beta - ksi*sqrt(m_N) / 2 * beta0;
				
		// calculate the fiber stress
		s = N*(2.0*Wl / J);

		//cout << n0.x << " " << n0.y << " " << n0.z << endl;

		// add the contribution from shear
		if (mu != 0.0)
		{
			mat3ds BmI = pt.LeftCauchyGreen() - mat3dd(1);
			s += (N*BmI).sym()*(mu / J);
		}

	}
	else
	{
		s.zero();
	}
	
	return s;
}

//-----------------------------------------------------------------------------
tens4ds FEFiberEntropyChain::FiberTangent(FEMaterialPoint& mp, const vec3d& n0)
{
	double ksi = m_ksi(mp);
	double mu = m_mu(mp);
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// deformation gradient
	mat3d &F = pt.m_F;
	double J = pt.m_J;
	
	mat3ds C = pt.RightCauchyGreen();
	tens4ds c;
	
	// fiber direction in global coordinate system
	//vec3d n0 = GetFiberVector(mp);
	
	// Calculate In = n0*C*n0
	double In = n0*(C*n0);

	// only take fibers in tension into consideration
	const double eps = m_epsf * std::numeric_limits<double>::epsilon();

	if ((In - 1.0) > eps)
	{
		// get the global spatial fiber direction in current configuration
		vec3d nt = F*n0;
		
		// calculate the outer product of nt
		mat3ds N = dyad(nt);
		tens4ds NxN = dyad1s(N);
		
		// calculate strain energy 2nd derivative
		const double f_a[] = { 3.0, 9.0 / 5.0, 297.0 / 175.0, 1539.0 / 875.0, 126117.0 / 67373.0, 43733439.0 / 21896875.0, \
			231321177.0 / 109484375.0, 20495009043.0 / 9306171875.0, 1073585186448381.0 / 476522530859375.0, 4387445039583.0 / 1944989921875.0, \
			1000263375846831627.0 / 453346207767578125.0, 280865021365240713.0 / 133337119931640625.0, 148014872740758343473.0 / 75350125192138671875.0, 137372931237386537808993.0 / 76480377070020751953125.0, \
			41722474198742657618857192737.0 / 25674386102028896409912109375.0, 12348948373636682700768301125723.0 / 8344175483159391333221435546875.0, 5001286000741585238340074032091091.0 / 3590449627018291035442047119140625.0, \
			185364329915163811141785118512534489.0 / 132846636199676768311355743408203125.0, 6292216384025878939310787532157558451.0 / 4160197291516193533960877227783203125.0, 299869254759556271677902570230858640837.0 / 170568088952163934892395966339111328125.0, \
			316689568216860631885141475537178451746044283.0 / 148810301691076651811854156590061187744140625.0, 670194310437429598283653289122392937145371137.0 / 258140319260030926612400067554187774658203125.0, 19697015384373759058671314622426656486031717197919.0 / 6332687088594936951180328352888571262359619140625.0, \
			178793788985653424246012689916144867915861856840849.0 / 49756827124674504616416865629838774204254150390625.0, 323844166067349737493036492206152479344269351967043143667.0 / 82039106319990447886052744467343800746748447418212890625.0, 200808116689754604893460969866238617668631975356485302537199.0 / 49409916306357883385918130190559334540655314922332763671875.0, \
			27506481209689719715275759452624078040221544551995885750037973.0 / 7164437864421893090958128877631103508395020663738250732421875.0, 16356939619211770477227805130221533318985185730316048281126247721.0 / 5122573073061653560035062147506239008502439774572849273681640625.0, 30976199222209837888906735596203249520053107250807132769871859115868101.0 / 14801698741072505317694781203956860833766632412399612367153167724609375.0, \
			519588001407316958447129785511020819131555326399179970047767492196701159.0 / 902903623205422824379381653441368510859764577156376354396343231201171875.0 };
		double alpha = sqrt(In) / sqrt(m_N);
		double alpha0 = 1 / sqrt(m_N);
		double alphaI = 1.0 / sqrt(In) / (2.0 * sqrt(m_N)); 
		double alphaII = -pow(In, (-3.0 / 2.0)) / (4.0 * sqrt(m_N));

		/////////////////////////////////////////////////////////////////////////////////
		double beta = 0, beta0 = 0, beta_alpha = 0, alpha2 = alpha * alpha, alpha02 = alpha0 * alpha0;
		int num_fa = sizeof(f_a) / sizeof(f_a[0]);
		if (m_term > num_fa) {
			m_term = num_fa;
		}
		beta = 1 + f_a[m_term - 1] / f_a[m_term - 2] * alpha2;
		beta0 = 1 + f_a[m_term - 1] / f_a[m_term - 2] * alpha02;
		beta_alpha = 1 + f_a[m_term - 1] / f_a[m_term - 2] * alpha2 * (2 * m_term - 1) / (2 * m_term - 3);
		for (int i = 2; i < m_term; i++) {
			beta = 1 + f_a[m_term - i] / f_a[m_term - i - 1] * alpha2 * beta;
			beta0 = 1 + f_a[m_term - i] / f_a[m_term - i - 1] * alpha02 * beta0;
			beta_alpha = 1 + f_a[m_term - i] / f_a[m_term - i - 1] * alpha2 * beta_alpha * (2 * m_term - 2 * i + 1) / (2 * m_term - 2 * i - 1);
		};
		beta = f_a[0] * alpha * beta;
		beta0 = f_a[0] * alpha0 * beta0;
		beta_alpha = f_a[0] * beta_alpha;
		//////////////////////////////////////////////////////////////////////////////////////


		double Wll = ksi * m_N * (beta_alpha*alphaI*alphaI + beta*alphaII);;
		
		// calculate the fiber tangent
		c = NxN*(4.0*Wll / J);
        
        // add the contribution from shear
		if (mu != 0.0)
		{
			mat3ds B = pt.LeftCauchyGreen();
			c += dyad4s(N,B)*(mu/J);
		}
	}
	else
	{
		c.zero();
	}
	
	return c;
}

//-----------------------------------------------------------------------------
//! Strain energy density
double FEFiberEntropyChain::FiberStrainEnergyDensity(FEMaterialPoint& mp, const vec3d& n0)
{
	double ksi = m_ksi(mp);
	double mu = m_mu(mp);
    double sed = 0.0;

	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// loop over all integration points
	mat3ds C = pt.RightCauchyGreen();
    mat3ds C2 = C.sqr();
	
	// fiber direction in global coordinate system
	//vec3d n0 = GetFiberVector(mp);
	
	// Calculate In = n0*C*n0
	double In = n0*(C*n0);

	// only take fibers in tension into consideration
	const double eps = m_epsf * std::numeric_limits<double>::epsilon();
	if ((In - 1.0) > eps)
	{
		// calculate strain energy density
		double alpha = sqrt(In) / sqrt(m_N);
		double beta = alpha * (3.0 - alpha * alpha) / (1.0 - alpha * alpha) - 0.5 * pow(alpha, (10.0 / 3.0)) + 3.0 * pow(alpha, 5.0)*(alpha - 0.76)*(alpha - 1.0);
		double alpha0 = 1 / sqrt(m_N);
		double beta0 = alpha0 * (3.0 - alpha0 * alpha0) / (1.0 - alpha0 * alpha0) - 0.5 * pow(alpha0, (10.0 / 3.0)) + 3.0 * pow(alpha0, 5.0)*(alpha0 - 0.76)*(alpha0 - 1.0);
		double alpha00 = ksi* sqrt(m_N) / 2 * beta0 + ksi*m_N*log(beta0 / sinh(beta0));
		sed = ksi*m_N * (alpha*beta + log(beta / (sinh(beta)))) - ksi*sqrt(m_N)/2 * beta0 * In - alpha00;
		
		// add the contribution from shear
		sed += mu*(n0*(C2*n0)-2.0*(In-1.0)-1.0)/4.0;
	}
	else
	{
		sed = 0; 
	}

    return sed;
}
//-----------------------------------------------------------------------------
// FEElasticFiberEntropyChain
//-----------------------------------------------------------------------------

// define the material parameters
BEGIN_FECORE_CLASS(FEElasticFiberEntropyChain, FEElasticFiberMaterial)
    ADD_PARAMETER(m_fib.m_N   , FE_RANGE_GREATER_OR_EQUAL(0.0), "N");
    ADD_PARAMETER(m_fib.m_ksi, FE_RANGE_GREATER_OR_EQUAL(0.0), "ksi")->setUnits(UNIT_PRESSURE);
    ADD_PARAMETER(m_fib.m_mu , FE_RANGE_GREATER_OR_EQUAL(0.0), "mu" )->setUnits(UNIT_PRESSURE);
    ADD_PARAMETER(m_fib.m_term, FE_RANGE_GREATER_OR_EQUAL(3), "n_term");
END_FECORE_CLASS();

