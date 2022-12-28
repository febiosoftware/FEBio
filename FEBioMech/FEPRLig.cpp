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
#include "FEPRLig.h"

// Define the material parameters
BEGIN_FECORE_CLASS(FEPRLig, FEElasticMaterial)
	ADD_PARAMETER(m_c1, "c1");
	ADD_PARAMETER(m_c2, "c2");
	ADD_PARAMETER(m_u,  "mu");
	ADD_PARAMETER(m_v0, "v0");
	ADD_PARAMETER(m_m,  "m");
	ADD_PARAMETER(m_k, "k");

	ADD_PROPERTY(m_Q, "mat_axis")->SetFlags(FEProperty::Optional);

END_FECORE_CLASS();

extern tens4ds material_to_spatial(tens4ds& C, mat3d& F);

//////////////////////////////////////////////////////////////////////
//								PRLig						
//////////////////////////////////////////////////////////////////////

FEPRLig::FEPRLig(FEModel* pfem) : FEElasticMaterial(pfem)
{
}

							///////////////////////////
//////////////////////////////Cauchy Stress Calculation ////////////////////////////////////////////////
							////////////////////////////

mat3ds FEPRLig::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	
	// obtain the deformation tensor F
	mat3d &F = pt.m_F;

	// calculate left and right Cauchy-Green tensors and their squares
	mat3ds b  =  pt.LeftCauchyGreen();
	mat3ds c  =  pt.RightCauchyGreen();				
	mat3ds b2 =  b.sqr();
	mat3ds c2 =  c.sqr();

	// Define the 2nd order identity tensor 
	mat3dd I(1);

	// get the local coordinate systems
	mat3d Q = GetLocalCS(mp);

	// declare initial material direction vector
	vec3d a0 = Q.col(0);

	// calculate the current material axis lam*a = F*a0;
	vec3d a = F*a0;

	// normalize material axis and store fiber stretch ; 
	double lam = a.unit();



//////////////////////////Strain invariants and Derivatives///////////////////////////////////////////////
	
	// define scalar strain invariants i1:15
	double i1 =  c.tr();
	double i2 =  0.5*(i1*i1 - c2.tr());
	double i3 =  c.det(); 
	double i4 =	 a0*(c*a0);
	double i5 =  a0*(c2*a0);

	
	// define common terms in the strain energy derivatives
	double rooti4 = sqrt(i4);
	double CT1 = i2 - i1*i4 + i5;
	double CT2 = -2*(m_m-m_v0);
	double CT3 = -1 + rooti4;
	double CT4 = exp(4*CT3*m_m);
	double CT4inv= 1/CT4;
	double CT5= pow(i4,-2*(m_m - m_v0));
	double CT5inv= 1/CT5;
	double CT6 = m_k*log(CT4*CT5*CT1)/CT1;
	double CTe = exp(m_c2*CT3*CT3);
	double CT8= m_k*(-CT4*i1*CT5 + 2*CT4*(1/rooti4)*CT5*CT1*m_m - 2*CT4*(1/i4)*CT5*CT1*(m_m - m_v0));
	

	// calculate the strain energy first derivatives with respect to each of the scalar invariants
	double w1 =  ((m_u)/2) - (i4*CT6); 

	double w2 =  CT6;

	double w3 = (-2*m_u)/(4*i3);
	
	double w4 =  m_c1*CTe*(CT3)/(2*rooti4) + (1/CT1)*CT4inv*CT5inv*CT8*CT6*CT1*(1/m_k);

	double w5 =  CT6;

	
	// calculate dyadic tensor terms found in stress calculation 

	// calculate a_i*a_j
	mat3ds aa = dyad(a);
	// calculate a_i*b_jk*a_k + b_ik*a_k*aj
	vec3d ba = b*a;
	mat3ds aobas = dyads(a, ba);
	// calculate J = sqrt(i3)
	double J = sqrt(i3);

	// calculate cauchy stress
	mat3ds s = 2/J*(i3*w3*I + (w1 + i1*w2)*b - w2*b2 + i4*w4*aa + i4*w5*aobas);
	
//  Stopping point debugging;  comment out when not testing
//	double WPR = 0;  

	// return stress
	return s;
}			



									//////////////////////////////
////////////////////////////////////Elasticity Tensor Calculation//////////////////////////////////////////////////////////////////////////
									///////////////////////////////

tens4ds FEPRLig::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// obtain the deformation tensor F
	mat3d &F = pt.m_F;


	// calculate left and right Cauchy-Green tensors, and their squares
	mat3ds b  =  pt.LeftCauchyGreen();
	mat3ds c  =  pt.RightCauchyGreen();				
	mat3ds b2 =  b.sqr();
	mat3ds c2 =  c.sqr();

	// define the 2nd order identity tensor 
	mat3dd I(1);

	// get the local coordinate systems
	mat3d Q = GetLocalCS(mp);

	// declare initial material direction vectors
	vec3d a0 = Q.col(0);

	// calculate the current material axis lam*a = F*a0;
	vec3d a = F*a0;

	// normalize material axis and store fiber stretch ; 
	double lam = a.unit();

	
	//////////////////////////////////Strain invariants and Derivatives///////////////////////////////////////

	// define scalar strain invariants i1:i5
	double i1 =  c.tr();
	double i2 =  0.5*(i1*i1 - c2.tr());
	double i3 =  c.det();
	double i4 =	 a0*(c*a0);
	double i5 =  a0*(c2*a0);

	// define common terms in the strain energy derivatives
	double rooti4 = sqrt(i4);
	double CT1 = i2 - i1*i4 + i5;
	double CT1_2= CT1*CT1;
	double CT2 = -2*(m_m-m_v0);
	double CT3 = -1 + rooti4;
	double CT4 = exp(4*CT3*m_m);
	double CT4inv= 1/CT4;
	double CT5= pow(i4,-2*(m_m - m_v0));
	double CT5inv= 1/CT5;
	double CT6 = m_k*log(CT4*CT5*CT1)/CT1;
	double CTe = exp(m_c2*CT3*CT3);
	double i4_2= i4*i4;
	double CT7= CT6*1/CT1;
	double CT8= m_k*(-CT4*i1*CT5 + 2*CT4*(1/rooti4)*CT5*CT1*m_m - 2*CT4*(1/i4)*CT5*CT1*(m_m - m_v0));
	double CT9 = exp(m_c2*CT3*CT3);

	// calculate the strain energy first derivatives with respect to each of the scalar invariants
	double w1 =  ((m_u)/2) - (i4*CT6); 

	double w2 =  CT6;

	double w3 = (-2*m_u)/(4*i3);
	
	double w4 =  m_c1*CTe*(CT3)/(2*rooti4) + (1/CT1)*CT4inv*CT5inv*CT8*CT6*CT1*(1/m_k);

	double w5 =  CT6;

	// calculate the second strain strain energy derivatives with respect to each of the scalar invariants
	double w11 = i4_2*m_k/CT1_2 - i4_2*CT7;

	double w12 = -i4*m_k/CT1_2 + i4*CT7;

	double w13 = 0;

	double w14 = -(1/CT1_2)*(1/CT4)*i4*(1/CT5)*CT8- i1*i4*CT7 - CT6;

	double w15 = w12;

	double w22 = m_k/CT1_2 - CT7;

	double w23 = 0;
	
	double w24 = 1/CT1_2*CT4inv*CT5inv*CT8 + i1*CT7;

	double w25 = w22;

	double w33 = (2*m_u)/(4*i3*i3);

	double w34 = 0;

	double w35 = 0;

	double w44 = (1/(4*i4_2*CT1_2))*(m_c1*CT9*rooti4*(1 + 2*m_c2*(rooti4 - 2*i4 + rooti4*rooti4*rooti4))*CT1_2 + 4*m_k*pow((-2*(i2 + i5)
					*(CT3*m_m + m_v0)+i1*i4*(1 + 2*CT3*m_m + 2*m_v0)),2)- 4*m_k*(-2*i1*i4*(i2 + i5)*((-2 + rooti4)*m_m + 2*m_v0)
					+ (i2 + i5)*(i2 + i5)*((-2+ rooti4)*m_m + 2*m_v0) + i1*i1*i4_2*(1 + (-2+ rooti4)*m_m + 2*m_v0))
					*log(CT4*CT5*CT1));
	 

	double w45 = 1/(CT1_2)*CT4inv*CT5inv*CT8 + 1/CT1*CT4inv*CT5inv*m_k*(2*CT4*(1/rooti4)*CT5*m_m - 2*CT4*(1/i4)*CT5*(m_m - m_v0))
				*CT6*CT1*(1/m_k) - 1/(CT1_2)*CT4inv*CT5inv*CT8*CT6*CT1*(1/m_k);
	
	double w55 = w22;

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	// calculate the second order dyadic terms in the elasticity tensor calculation
	mat3ds  a02   = dyad(a0); 
	vec3d	ctimesa0 = c*a0;
	mat3ds  di5dc = dyads(a0, ctimesa0);

	mat3ds cinv=c.inverse();

	mat3ds  aa        = dyad(a);
	vec3d   btimesa   = b*a;
	mat3ds  pf_didc   = dyads(btimesa,a);

	// calculate the fourth order dyadic terms in the elasticity tensor calculation
	tens4ds bb	      = dyad1s(b);
	tens4ds bb2       = dyad1s(b, b2);
	tens4ds b2b2      = dyad1s(b2);
	tens4ds pfi4      = dyad4s(b);    
	tens4ds bI        = dyad1s(b,I);
	tens4ds b2I       = dyad1s(b2,I);
	tens4ds II        = dyad1s(I);  
	tens4ds pf_dcindc = dyad4s(I); 
	tens4ds ba        = dyad1s(b,aa);
	tens4ds b2a       = dyad1s(b2,aa);
	tens4ds Ia		  = dyad1s(I,aa);
	tens4ds b_pfdi	  = dyad1s(pf_didc,b);
	tens4ds b2_pfdi   = dyad1s(pf_didc, b2);
	tens4ds I_pfdi    = dyad1s(pf_didc, I);
	tens4ds aa_pfdi   = dyad1s(pf_didc,aa);
	tens4ds a4        = dyad1s(aa);
	tens4ds pfdi_pfdi = dyad1s(pf_didc);
	tens4ds pf_ddi5dc2= dyad4s(aa,b); 


// Calculate Spatial Elasticity tensor
  tens4ds D = 4/sqrt(i3)*((w11 +2*i1*w12 + w2 + i1*i1*w22)*bb
					-(w12+i1*w22)*(bb2)
					+ w22*b2b2 - w2*pfi4
					+ (i3*w13 + i1*i3*w23)*(bI)
					- i3*w23*(b2I)
					+ (i3*w3 + i3*i3*w33)*(II)
					- w3*i3*(pf_dcindc)
					+ (i4*w14+i1*i4*w24)*(ba)
					- i4*w24*(b2a)
					+ i3*i4*w34*(Ia)
					+ (i4*w15 + i1*i4*w25)*(b_pfdi)
					- i4*w25*(b2_pfdi)
					+ i3*i4*w35*(I_pfdi)
					+ w45*i4*i4*(aa_pfdi)
					+ w44*i4*i4*(a4)
					+ i4*i4*w55*(pfdi_pfdi) + i4*w5*(pf_ddi5dc2)
					); 
					


  //Stopping point for stepping thru and checking stress calculation;  comment out when not testing
	// double PPP=1;


/////////////////alternate form of elasticity tensor in C^2 basis instead of Cinverse///////////////////////
/* 
mat3ds b3= b2*b;

tens4ds bb3		=dyad1s(b,b3);
tens4ds b2b3	=dyad1s(b2,b3);
tens4ds b3_pfdi	=dyad1s(pf_didc,b3);
tens4ds b3a		=dyad1s(b3,aa);
tens4ds b3b3 = dyad1s(b3);
tens4ds pfdc2dc = dyad4s(b2,b); 

tens4ds Dc2 =  (4/sqrt(i3))*((w11 +w2 + 2*i1*w12 + 2*i2*w13 + 2*i1*i2*w23 + i1*i1*w22 + i2*i2*w33 + i1*w3)*bb
	-(w3 + w12 + i1*w13 + i2*w23 + i1*w22 + i1*i2*w33 + i1*i1*w23)*bb2
	+(w22 + 2*i1*w23 + i1*i1*w33)*b2b2
	- (w2 + i1*w3)*pfi4
	+ (w13 + i1*w23 + i2*w33)*bb3
	- (w23+i1*w33)*b2b3
	+ w33*b3b3 + w3*pfdc2dc
	+ (i4*w14 + i1*i4*w24 + i2*i4*w34)*ba
	- (i4*w24 + i1*i4*w34)*b2a
	+ i4*w34*b3a
	+ (i4*w15 + i1*i4*w25 + i2*i4*w35)*b_pfdi
	- (i4*w25 + i1*i4*w35)*b2_pfdi
	+ i4*i4*w44*a4
	+ i4*w35*b3_pfdi
	+ i4*i4*w45*aa_pfdi
	+ i4*i4*w55*(pfdi_pfdi) + i4*w5*(pf_ddi5dc2));


 double PPP2=0;
 */
	return D;
}

//! calculate strain energy density at material point
double FEPRLig::StrainEnergyDensity(FEMaterialPoint& pt)
{
    double sed = 0;

    FEElasticMaterialPoint& mp = *pt.ExtractData<FEElasticMaterialPoint>();
    
    // obtain the deformation tensor F
    mat3d &F = mp.m_F;
    
    // calculate left and right Cauchy-Green tensors and their squares
    mat3ds c  =  mp.RightCauchyGreen();
    mat3ds c2 =  c.sqr();
    
    // Define the 2nd order identity tensor
    mat3dd I(1);
    
    // get the local coordinate systems
    mat3d Q = GetLocalCS(pt);
    
    // declare initial material direction vector
    vec3d a0 = Q.col(0);
    
    // calculate the current material axis lam*a = F*a0;
    vec3d a = F*a0;
    
    // normalize material axis and store fiber stretch ;
    double lam = a.unit();
    
    //////////////////////////Strain invariants ///////////////////////////////////////////////
    
    // define scalar strain invariants i1:15
    double i1 =  c.tr();
    double i2 =  0.5*(i1*i1 - c2.tr());
    double i3 =  c.det();
    double i4 =     a0*(c*a0);
    double i5 =  a0*(c2*a0);

    double Wfiber = 0.5*m_c1/m_c2*(exp(m_c2*pow(lam-1,2))-1);
    double Wmatrix = m_u/2*(i1-3) - m_u*log(sqrt(i3));
    double Wvol = m_k/2*pow(log((i5-i1*i4+i2)/(pow(i4, 2*(m_m-m_v0)*exp(-4*m_m*(lam-1))))),2);
    
    sed = Wfiber + Wmatrix + Wvol;
    
    return sed;
}
