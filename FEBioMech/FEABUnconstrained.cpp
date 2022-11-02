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


#include "FECore/stdafx.h"
#include "FEABUnconstrained.h"
#include <iostream>
#include <typeinfo>
using namespace std;

BEGIN_FECORE_CLASS(FEABUnconstrained, FEElasticMaterial);
	ADD_PARAMETER(m_ksi, FE_RANGE_GREATER(0.0), "ksi");
	ADD_PARAMETER(m_N, FE_RANGE_GREATER(0.0), "N");
	ADD_PARAMETER(m_term, FE_RANGE_CLOSED(3,30), "n_term");
	ADD_PARAMETER(m_kappa, FE_RANGE_GREATER(0.0), "kappa");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FEABUnconstrained::FEABUnconstrained(FEModel* pfem) : FEElasticMaterial(pfem)
{
	m_N = 100.0; 
	m_ksi = 0.00001;
	m_kappa = 0.0;
	m_term = 30;
	m_eps = 1e-12;
}

//-----------------------------------------------------------------------------
void FEABUnconstrained::EigenValues(mat3ds& A, double l[3], vec3d r[3], const double eps)
{
	A.eigen(l, r);
	
	// correct for numerical inaccuracy
	double d01 = fabs(l[0] - l[1]);
	double d12 = fabs(l[1] - l[2]);
	double d02 = fabs(l[0] - l[2]);
	
	if (d01 < eps) l[1] = l[0]; //= 0.5*(l[0]+l[1]);
	if (d02 < eps) l[2] = l[0]; //= 0.5*(l[0]+l[2]);
	if (d12 < eps) l[2] = l[1]; //= 0.5*(l[1]+l[2]);
	
}

//-----------------------------------------------------------------------------
//! Calculates the Cauchy stress
mat3ds FEABUnconstrained::Stress(FEMaterialPoint& mp)
{
	double ksi = m_ksi(mp);
	double kappa = m_kappa(mp);
	// extract elastic material data
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// jacobian
	double J = pt.m_J;

	// get the left Cauchy-Green tensor
	mat3ds b = pt.LeftCauchyGreen();

	// get the eigenvalues and eigenvectors of b
	double lam2[3];	// these are the squares of the eigenvalues of V
	vec3d ev[3];
	EigenValues(b, lam2, ev, m_eps);

	// get the eigenvalues of V
	double lam[3];
	mat3ds N[3];
	for (int i = 0; i < 3; ++i) {
		lam[i] = sqrt(lam2[i]);
		N[i] = dyad(ev[i]);
	}

	// Invariants of b
	double I1 = b.tr();

	// stress
	mat3ds s;
	s.zero();

	const double f_a[] = { 3.0, 9.0 / 5.0, 297.0 / 175.0, 1539.0 / 875.0, 126117.0 / 67373.0, 43733439.0 / 21896875.0, \
		231321177.0 / 109484375.0, 20495009043.0 / 9306171875.0, 1073585186448381.0 / 476522530859375.0, 4387445039583.0 / 1944989921875.0, \
			1000263375846831627.0 / 453346207767578125.0, 280865021365240713.0 / 133337119931640625.0, 148014872740758343473.0 / 75350125192138671875.0, 137372931237386537808993.0 / 76480377070020751953125.0, \
			41722474198742657618857192737.0 / 25674386102028896409912109375.0, 12348948373636682700768301125723.0 / 8344175483159391333221435546875.0, 5001286000741585238340074032091091.0 / 3590449627018291035442047119140625.0, \
			185364329915163811141785118512534489.0 / 132846636199676768311355743408203125.0, 6292216384025878939310787532157558451.0 / 4160197291516193533960877227783203125.0, 299869254759556271677902570230858640837.0 / 170568088952163934892395966339111328125.0, \
			316689568216860631885141475537178451746044283.0 / 148810301691076651811854156590061187744140625.0, 670194310437429598283653289122392937145371137.0 / 258140319260030926612400067554187774658203125.0, 19697015384373759058671314622426656486031717197919.0 / 6332687088594936951180328352888571262359619140625.0, \
			178793788985653424246012689916144867915861856840849.0 / 49756827124674504616416865629838774204254150390625.0, 323844166067349737493036492206152479344269351967043143667.0 / 82039106319990447886052744467343800746748447418212890625.0, 200808116689754604893460969866238617668631975356485302537199.0 / 49409916306357883385918130190559334540655314922332763671875.0, \
			27506481209689719715275759452624078040221544551995885750037973.0 / 7164437864421893090958128877631103508395020663738250732421875.0, 16356939619211770477227805130221533318985185730316048281126247721.0 / 5122573073061653560035062147506239008502439774572849273681640625.0, 30976199222209837888906735596203249520053107250807132769871859115868101.0 / 14801698741072505317694781203956860833766632412399612367153167724609375.0, \
			519588001407316958447129785511020819131555326399179970047767492196701159.0 / 902903623205422824379381653441368510859764577156376354396343231201171875.0 };
	double alpha = sqrt(I1) / sqrt(3.0 * m_N);
	double alpha0 = 1.0 / sqrt(m_N);

	///////////////////////////////
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
	};
	beta = f_a[0] * alpha * beta;
	beta0 = f_a[0] * alpha0 * beta0;
	////////////////////////////////

	double WI = ksi * m_N * (beta / sqrt(I1) / 2.0 / sqrt(3.0 * m_N));
	double W1[3];
	double T[3], T_NH[3];
	double I1_1[3];
	double W10 = ksi * m_N * (beta0 / sqrt(3.0) / 2.0 / sqrt(3.0 * m_N)) * 2.0;

	/////////////////////////////////////////////////////////////////////
	double aaaI = exp((I1) / 3 - 1); //Simple

	//double TCom = kappa * log(J) / J; //NH Term
	double TCom = kappa * (J * J - 1) / J; //NH Term

	for (int i = 0; i < 3; ++i) {
		I1_1[i] = 2 * lam[i];
		W1[i] = WI * I1_1[i] - W10 / lam[i];
		T[i] = lam[i] * W1[i] / J;
		T[i] += TCom; //NH Term 
		s += N[i] * T[i];
	}


	return s;
}

//-----------------------------------------------------------------------------
//! Calculates the spatial tangent
tens4ds FEABUnconstrained::Tangent(FEMaterialPoint& mp)
{
	double ksi = m_ksi(mp);
	double kappa = m_kappa(mp);
	int i, j;

	// extract elastic material data
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// jacobian
	double J = pt.m_J;

	// get the left Cauchy-Green tensor
	mat3ds b = pt.LeftCauchyGreen();

	// Invariants of b
	double I1 = b.tr();

	// get the eigenvalues and eigenvectors of b
	double lam2[3];	// these are the squares of the eigenvalues of V
	vec3d ev[3];
	EigenValues(b, lam2, ev, m_eps);

	// get the eigenvalues of V
	double lam[3];
	mat3ds N[3];
	for (i = 0; i < 3; ++i) {
		lam[i] = sqrt(lam2[i]);
		N[i] = dyad(ev[i]);
	}

	// principal stresses
	mat3ds s;
	s.zero();

	const double f_a[] = { 3.0, 9.0 / 5.0, 297.0 / 175.0, 1539.0 / 875.0, 126117.0 / 67373.0, 43733439.0 / 21896875.0, \
		231321177.0 / 109484375.0, 20495009043.0 / 9306171875.0, 1073585186448381.0 / 476522530859375.0, 4387445039583.0 / 1944989921875.0, \
			1000263375846831627.0 / 453346207767578125.0, 280865021365240713.0 / 133337119931640625.0, 148014872740758343473.0 / 75350125192138671875.0, 137372931237386537808993.0 / 76480377070020751953125.0, \
			41722474198742657618857192737.0 / 25674386102028896409912109375.0, 12348948373636682700768301125723.0 / 8344175483159391333221435546875.0, 5001286000741585238340074032091091.0 / 3590449627018291035442047119140625.0, \
			185364329915163811141785118512534489.0 / 132846636199676768311355743408203125.0, 6292216384025878939310787532157558451.0 / 4160197291516193533960877227783203125.0, 299869254759556271677902570230858640837.0 / 170568088952163934892395966339111328125.0, \
			316689568216860631885141475537178451746044283.0 / 148810301691076651811854156590061187744140625.0, 670194310437429598283653289122392937145371137.0 / 258140319260030926612400067554187774658203125.0, 19697015384373759058671314622426656486031717197919.0 / 6332687088594936951180328352888571262359619140625.0, \
			178793788985653424246012689916144867915861856840849.0 / 49756827124674504616416865629838774204254150390625.0, 323844166067349737493036492206152479344269351967043143667.0 / 82039106319990447886052744467343800746748447418212890625.0, 200808116689754604893460969866238617668631975356485302537199.0 / 49409916306357883385918130190559334540655314922332763671875.0, \
			27506481209689719715275759452624078040221544551995885750037973.0 / 7164437864421893090958128877631103508395020663738250732421875.0, 16356939619211770477227805130221533318985185730316048281126247721.0 / 5122573073061653560035062147506239008502439774572849273681640625.0, 30976199222209837888906735596203249520053107250807132769871859115868101.0 / 14801698741072505317694781203956860833766632412399612367153167724609375.0, \
			519588001407316958447129785511020819131555326399179970047767492196701159.0 / 902903623205422824379381653441368510859764577156376354396343231201171875.0 };
	double alpha = sqrt(I1) / sqrt(3.0 * m_N);
	double alpha0 = 1.0 / sqrt(m_N);

		//////////////////////////////////////
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
	////////////////////////////////////////

	double WI = ksi * m_N * (beta / sqrt(I1) / 2.0 / sqrt(3.0 * m_N));
	double WII = ksi * m_N * (beta_alpha / I1 / sqrt(3.0 * m_N) - beta * pow(I1, (-3.0 / 2.0))) / 4.0 / sqrt(3.0 * m_N);
	double W1[3];
	double W11[3][3];
	double T[3], T_NH[3];
	double I1_1[3];
	double I1_11 = 2;
	double W10 = ksi * m_N * (beta0 / sqrt(3.0) / 2.0 / sqrt(3 * m_N)) * 2.0;

	/////////////////////////////////////////////////////////////////////

	double TCom = kappa * (J * J - 1) / J; //NH Term

	for (int i = 0; i < 3; ++i) {
		I1_1[i] = 2 * lam[i];
		W1[i] = WI * I1_1[i] - W10 / lam[i];
		T[i] = lam[i] * W1[i] / J;
		T[i] += TCom;
		s += N[i] * T[i];
	}

	// coefficients appearing in elasticity tensor
	double D[3][3], E[3][3];

	double D_NH[3][3], E_NH[3][3]; //NH Term     

	double DComii = 2 * kappa / J + TCom; //NH Term
	double DComij = 2 * kappa * J; //NH Term
	double EComij = 2 * (-kappa * J * J + kappa) / J; //NH Term

	for (i = 0; i < 3; ++i) {

		I1_1[i] = 2 * lam[i];
		W11[i][i] = WII * I1_1[i] * I1_1[i] + WI * I1_11 + W10 / lam2[i];
		D[i][i] = lam2[i] * W11[i][i] / J - T[i];
		D[i][i] += DComii;
		for (j = i + 1; j < 3; ++j) {

			W11[i][j] = WII * I1_1[i] * I1_1[j];
			D[i][j] = lam[i] * lam[j] * W11[i][j] / J;
			D[i][j] += DComij;

			if (lam2[j] != lam2[i]) {

				E[i][j] = 2 * (lam2[j] * T[i] - lam2[i] * T[j]) / (lam2[i] - lam2[j]);
			}
			else {
				E[i][j] = 2 * W10 / J;
				E[i][j] += EComij;
			}
		}
	}

	// spatial elasticity tensor
	tens4ds c(0.0);
	mat3dd I(1.0);
	for (i = 0; i < 3; ++i) {
		c += dyad1s(N[i]) * D[i][i];
		for (j = i + 1; j < 3; ++j) {
			c += dyad1s(N[i], N[j]) * D[i][j];
			c += dyad4s(N[i], N[j]) * E[i][j];
		}
	}

	return c;
}

//-----------------------------------------------------------------------------
double FEABUnconstrained::StrainEnergyDensity(FEMaterialPoint& mp)
{
	double ksi = m_ksi(mp);
	double kappa = m_kappa(mp);
	// extract elastic material data
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// jacobian
	double J = pt.m_J;
    double lnJ = log(J);
	
	// get the left Cauchy-Green tensor
	mat3ds b = pt.LeftCauchyGreen();
	
	// get the eigenvalues and eigenvectors of b
	double lam2[3];	// these are the squares of the eigenvalues of V
	vec3d ev[3];
	EigenValues(b, lam2, ev, m_eps);
	
	// get the eigenvalues of V
	double lam[3];
	lam[0] = sqrt(lam2[0]);
	lam[1] = sqrt(lam2[1]);
	lam[2] = sqrt(lam2[2]);
	
	// Invariants of b
	double I1 = b.tr();

	// strain energy density

	double alpha = sqrt(I1) / sqrt(3*m_N);
	double beta = alpha * (3.0 - alpha * alpha) / (1.0 - alpha * alpha) -0.5 * pow(alpha, (10.0 / 3.0)) + 3.0 * pow(alpha, 5.0)*(alpha - 0.76)*(alpha - 1.0);
	double alpha0 = 1 / sqrt(m_N);
	double beta0 = alpha0 * (3.0 - alpha0 * alpha0) / (1.0 - alpha0 * alpha0) -0.5 * pow(alpha0, (10.0 / 3.0)) + 3.0 * pow(alpha0, 5.0)*(alpha0 - 0.76)*(alpha0 - 1.0);
	
	double W10 = ksi * m_N * (beta0 / sqrt(3) / 2 / sqrt(3 * m_N)) * 2;
	double sedCom = kappa / kappa / kappa * cosh(kappa*(J - 1) - 1);

	double sed = ksi * m_N * (sqrt(I1)*beta / sqrt(3 * m_N) + log(beta / sinh(beta))) - W10 * lnJ +sedCom;
	
	return sed;
}
