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
#include "FEMRVonMisesFibers.h"
#include <fstream>

// define the material parameters
BEGIN_FECORE_CLASS(FEMRVonMisesFibers, FEUncoupledMaterial)
	ADD_PARAMETER(c1, "c1");
	ADD_PARAMETER(c2, "c2");
	ADD_PARAMETER(m_fib.m_c3, "c3");
	ADD_PARAMETER(m_fib.m_c4, "c4");
	ADD_PARAMETER(m_fib.m_c5, "c5");
	ADD_PARAMETER(m_fib.m_lam1, "lam_max");
	// Fiber Concentration Factor
	ADD_PARAMETER(kf, "kf");
	// Preferred fiber orientation IN RADIANS
	ADD_PARAMETER(tp, "tp");
	// Number of Gauss Integration Points 
	ADD_PARAMETER(gipt, "gipt");
	// Choice of von Mises distribution; 1: semi-circular von Mises distribution; = 2: constrained von Mises distribution
	ADD_PARAMETER(vmc, "vmc"); 
	// Exponent for the constrained von Mises distribution
	ADD_PARAMETER(var_n, "var_n"); 

	ADD_PROPERTY(m_Q, "mat_axis")->SetFlags(FEProperty::Optional);

END_FECORE_CLASS();

//=============================================================================
FEMRVonMisesMaterialPoint::FEMRVonMisesMaterialPoint(FEMaterialPointData* mp) : FEMaterialPointData(mp)
{

}

FEMaterialPointData* FEMRVonMisesMaterialPoint::Copy()
{
	FEMRVonMisesMaterialPoint* pt = new FEMRVonMisesMaterialPoint(*this);
	if (m_pNext) pt->m_pNext = m_pNext->Copy();
	return pt;
}

//-----------------------------------------------------------------------------
void FEMRVonMisesMaterialPoint::Serialize(DumpStream& ar)
{
	FEMaterialPointData::Serialize(ar);
	ar & m_kf & m_tp;
}

//////////////////////////////////////////////////////////////////////
// FEVonMisesFibers
/////////////////////////////////////////////////////////////////////
//
// This material model has been implemented by Cecile Gouget and Michael Girard
// 
// This is a fiber-based model intended to describe the behavior of the sclera - the white coat of the eye - and other thin soft tissues. 
// Instead of having a single fiber orientation (Transverse Isotropy), we consider that collagen fiber alignment is 
// multi-directional at local material points, confined within the plane tangent to the tissue surface, and described by either: 
// 
// 1) A semi-circular von Mises distribution function. Such model was originally implemented within Nike3d [1].
// 2) A constrained von Mises distribution function [2]
//
// This latter function is formed as a weighted mixture of the semi-circular uniform distribution and the semi-circular von Mises distribution. 
// Such function has been shown to better describe the planar anisotropy of thin soft tissues and was validated using small angle light scattering measurements 
// tissue anisotropy. It is therefore the preferred function to use. 
//
// REFERENCES
// ------------
// [1] Girard MJA, Downs JC, Burgoyne CF, and Suh JKF, 2009, "Peripapillary and posterior scleral mechanics - Part I: 
// Development of an anisotropic hyperelastic constitutive model," ASME J. Biomech. Eng., 131, p. 051011.  
// ------------
// [2] Gouget CLM, Girard MJA, Ethier CR, 2012, "A constrained von Mises distribution to describe fiber organization 
// in thin soft tissues," Biomechanics And Modeling in Mechanobiolgy, 11(3-4), p. 475-482. 
//
//
/////////////////////////////////////////////////////////////////////

// This is the modified Bessel function of the first kind (order zero)
// It is used to compute the probability distribution function 
// of the fibers - the semi-circular von-Mises distribution 
double bessi0(double X)
{
	double Y,P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,AX,BX;
	P1=1.0; P2=3.5156229; P3=3.0899424; P4=1.2067429;
	P5=0.2659732; P6=0.360768e-1; P7=0.45813e-2;
	Q1=0.39894228; Q2=0.1328592e-1; Q3=0.225319e-2;
	Q4=-0.157565e-2; Q5=0.916281e-2; Q6=-0.2057706e-1;
	Q7=0.2635537e-1; Q8=-0.1647633e-1; Q9=0.392377e-2;
	if (fabs(X) < 3.75)
	{
		Y=(X/3.75)*(X/3.75);
		return (P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))));
	}
	else
	{
		AX=fabs(X);
		Y=3.75/AX;
		BX=exp(AX)/sqrt(AX);
		AX=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))));
		return (AX*BX);
	}
}

// This is the modified Bessel function of the first kind (order one)
// It is used to compute the probability distribution function 
// of the fibers - the constrained von-Mises distribution 
double bessi1(double X) 
{
	double Y,P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,AX,BX;
	P1=0.5; P2=0.87890594; P3=0.51498869; P4=0.15084934;
	P5=0.2658733e-1; P6=0.301532e-2; P7=0.32411e-3;
	Q1=0.39894228; Q2=-0.3988024e-1; Q3=-0.362018e-2;
	Q4=0.163801e-2; Q5=-0.1031555e-1; Q6=0.2282967e-1;
	Q7=-0.2895312e-1; Q8=0.1787654e-1; Q9=-0.420059e-2;
	if (fabs(X) < 3.75) 
	{
        Y=(X/3.75)*(X/3.75);
        return(X*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))));
	}
	else 
	{
        AX=fabs(X);
        Y=3.75/AX;
        BX=exp(AX)/sqrt(AX);
        AX=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))));
        return (AX*BX);
	}
}

//-----------------------------------------------------------------------------
FEMRVonMisesFibers::FEMRVonMisesFibers(FEModel* pfem) : FEUncoupledMaterial(pfem), m_fib(pfem) 
{
	m_fib.SetParent(this);
}

//-----------------------------------------------------------------------------
FEMaterialPointData* FEMRVonMisesFibers::CreateMaterialPointData()
{
	FEMRVonMisesMaterialPoint* pt = new FEMRVonMisesMaterialPoint(FEUncoupledMaterial::CreateMaterialPointData());
	pt->m_kf = kf;
	pt->m_tp = tp;
	return pt;
}

//-----------------------------------------------------------------------------
mat3ds FEMRVonMisesFibers::DevStress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	FEMRVonMisesMaterialPoint& vmpt = *mp.ExtractData<FEMRVonMisesMaterialPoint>();

	// get some of the material parmaters from the material point
	double kf = vmpt.m_kf;
	double tp = vmpt.m_tp;

	// deformation gradient
	mat3d& F = pt.m_F;
	double J = pt.m_J;

	//// FIRST, WE CALCULATE THE MOONEY RIVLIN PART OF THE STRESS
	// calculate deviatoric left Cauchy-Green tensor
	mat3ds B = pt.DevLeftCauchyGreen();
	
	// calculate square of B
	mat3ds B2 = B.sqr();

	// Invariants of B (= invariants of C)
	// Note that these are the invariants of Btilde, not of B!
	double I1 = B.tr();
	double I2 = 0.5*(I1*I1 - B2.tr());

	// --- TODO: put strain energy derivatives here ---
	// Wi = dW/dIi
	double W1 = c1;
	double W2 = c2;

	// calculate T = F*dW/dC*Ft
	mat3ds T = B*(W1 + W2*I1) - B2*W2;
	
	// calculate stress s = 2/J * dev(T) 
	mat3ds s =  T.dev()*(2.0/J);

	//NOW, WE CALCULATE THE FIBER COMPONENT OF THE STRESS by the mean of a Gauss integration
	const double pi = 3.14159265358979323846;
	vec3d a0, b0, c0;			// current local material axis
	vector<double> pp(gipt), wmg(gipt);		//gipt Gauss points and weights
	double aa, bb, pas;
	int i;
	double distribution;			//experimental distribution function of fiber orientations
//	double sigma;
	double beta;
 
        //std::cout << gipt << "\n"; 
	
	//-=-=- Gauss Points -=-=- 
	// Compute angles for Gaussian Quadrature  
	// for a function f, a Gaussian quadrature rule follows 
	//
	//                         n
	//    /a                 -----
	//   |           b - a   \             b - a      a + b
	//   | f(t)dt = -------   \    wi * f( ----- xi + ----- ) 
	//   |             2      /              2          2
	//  /  b                 /
	//                       -----    
	//                       i = 1 
	//
	// where xi are the points and wi are the weights 
	//
	// For n = 2, we have 
	// w1 = w2 = 1.0
	// x1 = - x2 = 1/sqrt(3.0)
	//
	// In the following we used "gipt" Gauss points for the numerical integration
	//-=-=-=-=-=-=-=-=-=-=-=-=
       
	// step
	pas = pi / (gipt/10);

	//Gauss Points
	for(int j=0; j<gipt; j+=10)
	{
         aa = tp - pi/2.0 +  j/10    * pas;
         bb = tp - pi/2.0 + (j/10+1) * pas; 
         pp[j+0] = pas/2.0 * ( + 0.14887434 ) + (aa+bb)/2.0;
         pp[j+1] = pas/2.0 * ( - 0.14887434 ) + (aa+bb)/2.0;
         pp[j+2] = pas/2.0 * ( + 0.43339539 ) + (aa+bb)/2.0;
         pp[j+3] = pas/2.0 * ( - 0.43339539 ) + (aa+bb)/2.0;
         pp[j+4] = pas/2.0 * ( + 0.67940957 ) + (aa+bb)/2.0;
         pp[j+5] = pas/2.0 * ( - 0.67940957 ) + (aa+bb)/2.0;
         pp[j+6] = pas/2.0 * ( + 0.86506337 ) + (aa+bb)/2.0;
         pp[j+7] = pas/2.0 * ( - 0.86506337 ) + (aa+bb)/2.0;
         pp[j+8] = pas/2.0 * ( + 0.97390653 ) + (aa+bb)/2.0;
         pp[j+9] = pas/2.0 * ( - 0.97390653 ) + (aa+bb)/2.0;
 	}

	//Weights for Gaussian Quadrature  
	for(int j=0; j<gipt; j+=10)
	{
		wmg[j+0] = 0.29552422 * pas/2.0;        
         wmg[j+1] = 0.29552422 * pas/2.0;       
         wmg[j+2] = 0.26926672 * pas/2.0;
         wmg[j+3] = 0.26926672 * pas/2.0;
         wmg[j+4] = 0.21908636 * pas/2.0;
         wmg[j+5] = 0.21908636 * pas/2.0;
         wmg[j+6] = 0.14945135 * pas/2.0;
         wmg[j+7] = 0.14945135 * pas/2.0;
         wmg[j+8] = 0.06667134 * pas/2.0;
         wmg[j+9] = 0.06667134 * pas/2.0;
 	}

	// get the local coordinate systems
	mat3d Q = GetLocalCS(mp);

	// get the initial fiber direction : a0 is the fiber direction,  
	// b0 a vector of the plane of the fibers perpendicular to a0
	// (specified in mat_axis : a0 = unit(a) ; b0 = unit((a0^d0)^a0)
	a0.x = Q[0][0]; b0.x = Q[0][1]; c0.x = Q[0][2];
	a0.y = Q[1][0]; b0.y = Q[1][1]; c0.y = Q[1][2];
	a0.z = Q[2][0]; b0.z = Q[2][1]; c0.z = Q[2][2];

	// Think of it as #gipt different fiber orientations that are distributed within a plane
	// of coordinate system (a0,b0)
	for(i=0;i<gipt;++i)
	{
		// vector describing the fibers which make an angle pp[i] with vector a0 in plane (a0,b0)
		vec3d n0;
		n0.x = cos(pp[i])*a0.x+sin(pp[i])*b0.x;
		n0.y = cos(pp[i])*a0.y+sin(pp[i])*b0.y;
		n0.z = cos(pp[i])*a0.z+sin(pp[i])*b0.z;

		// probability of having a fiber along this vector 
		if (vmc==1) // Semi-circular von Mises distribution
		{
			distribution = 1. / ( pi * bessi0(kf) ) * exp( kf * cos( 2.* (pp[i]-tp) ) );
		}
		else if (vmc==2) // Constrained von Mises distribution
		{
			beta = pow(bessi1(kf)/bessi0(kf),var_n);  // Weighting factor (isotropic vs anisotropic)
			distribution = (1.-beta)/pi + beta / ( pi * bessi0(kf) ) * exp( kf * cos( 2.* (pp[i]-tp) ) );
		}		
		
		// add the fiber stress by Gauss integration : m_fib.Stress(mp) is the deviatoric stress due to the fibers in direction pp[i], 
		// distribution is the probability of having a fiber in this direction, and wmg[i] is the Gauss Weight for integration
		s += wmg[i]*distribution*m_fib.DevFiberStress(mp, n0);
	}

	return s;
}

tens4ds FEMRVonMisesFibers::DevTangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	FEMRVonMisesMaterialPoint& vmpt = *mp.ExtractData<FEMRVonMisesMaterialPoint>();

	// get some of the material parmaters from the material point
	double kf = vmpt.m_kf;
	double tp = vmpt.m_tp;

	// deformation gradient
	mat3d& F = pt.m_F;
	double J = pt.m_J;
	double Ji = 1.0/J;

	//// FIRST, WE CALCULATE THE MOONEY RIVLIN PART OF THE TANGENT
	// calculate deviatoric left Cauchy-Green tensor: B = F*Ft
	mat3ds B = pt.DevLeftCauchyGreen();

	// calculate square of B
	mat3ds B2 = B.sqr();

	// Invariants of B (= invariants of C)
	double I1 = B.tr();
	double I2 = 0.5*(I1*I1 - B2.tr());

	// --- TODO: put strain energy derivatives here ---
	// Wi = dW/dIi
	double W1, W2;
	W1 = c1;
	W2 = c2;

	// calculate dWdC:C
	double WC = W1*I1 + 2*W2*I2;

	// calculate C:d2WdCdC:C
	double CWWC = 2*I2*W2;

	mat3dd I(1);	// Identity
	tens4ds IxI = dyad1s(I);
	tens4ds I4  = dyad4s(I);
	tens4ds BxB = dyad1s(B);
	tens4ds B4  = dyad4s(B);

	// deviatoric cauchy-stress, trs = trace[s]/3
	mat3ds devs = pt.m_s.dev();

	// d2W/dCdC:C
	mat3ds WCCxC = B*(W2*I1) - B2*W2;

	tens4ds cw = (BxB - B4)*(W2*4.0*Ji) - dyad1s(WCCxC, I)*(4.0/3.0*Ji) + IxI*(4.0/9.0*Ji*CWWC);
	tens4ds c = dyad1s(devs, I)*(2.0/3.0) + (I4 - IxI/3.0)*(4.0/3.0*Ji*WC) + cw;

	//NOW, WE CALCULATE THE FIBERS CONTRIBUTION TO THE TANGENT by the mean of a Gauss integration
	// current local material axis
	const double pi = 3.14159265358979323846;	
	vec3d a0, b0, c0;
	vector<double> pp(gipt), wmg(gipt);
	double aa, bb, pas;
	int i;
	double distribution;
//	double sigma;
	double beta;
	
	//-=-=- Gauss Points -=-=- 
	// Compute angles for Gaussian Quadrature  
	// for a function f, a Gaussian quadrature rule follows 
	//
	//                         n
	//    /a                 -----
	//   |           b - a   \             b - a      a + b
	//   | f(t)dt = -------   \    wi * f( ----- xi + ----- ) 
	//   |             2      /              2          2
	//  /  b                 /
	//                       -----    
	//                       i = 1 
	//
	// where xi are the points and wi are the weights 
	//
	// For n = 2, we have 
	// w1 = w2 = 1.0
	// x1 = - x2 = 1/sqrt(3.0)
	//
	// In the following we used "gipt" Gauss points for the numerical integration
	//-=-=-=-=-=-=-=-=-=-=-=-=
       
	// step
	pas = pi / (gipt/10);

	//Gauss Points
	for(int j=0; j<gipt; j+=10)
	{
         aa = tp - pi/2.0 +  j/10    * pas;
         bb = tp - pi/2.0 + (j/10+1) * pas; 
         pp[j+0] = pas/2.0 * ( + 0.14887434 ) + (aa+bb)/2.0;
         pp[j+1] = pas/2.0 * ( - 0.14887434 ) + (aa+bb)/2.0;
         pp[j+2] = pas/2.0 * ( + 0.43339539 ) + (aa+bb)/2.0;
         pp[j+3] = pas/2.0 * ( - 0.43339539 ) + (aa+bb)/2.0;
         pp[j+4] = pas/2.0 * ( + 0.67940957 ) + (aa+bb)/2.0;
         pp[j+5] = pas/2.0 * ( - 0.67940957 ) + (aa+bb)/2.0;
         pp[j+6] = pas/2.0 * ( + 0.86506337 ) + (aa+bb)/2.0;
         pp[j+7] = pas/2.0 * ( - 0.86506337 ) + (aa+bb)/2.0;
         pp[j+8] = pas/2.0 * ( + 0.97390653 ) + (aa+bb)/2.0;
         pp[j+9] = pas/2.0 * ( - 0.97390653 ) + (aa+bb)/2.0;
 	}

	//Weights for Gaussian Quadrature  
	for(int j=0; j<gipt; j+=10)
	{
		wmg[j+0] = 0.29552422 * pas/2.0;        
		wmg[j+1] = 0.29552422 * pas/2.0;       
		wmg[j+2] = 0.26926672 * pas/2.0;
		wmg[j+3] = 0.26926672 * pas/2.0;
		wmg[j+4] = 0.21908636 * pas/2.0;
		wmg[j+5] = 0.21908636 * pas/2.0;
		wmg[j+6] = 0.14945135 * pas/2.0;
		wmg[j+7] = 0.14945135 * pas/2.0;
		wmg[j+8] = 0.06667134 * pas/2.0;
		wmg[j+9] = 0.06667134 * pas/2.0;
 	}

	// get the local coordinate systems
	mat3d Q = GetLocalCS(mp);

	// get the initial fiber direction
	a0.x = Q[0][0]; b0.x = Q[0][1]; c0.x = Q[0][2];
	a0.y = Q[1][0]; b0.y = Q[1][1]; c0.y = Q[1][2];
	a0.z = Q[2][0]; b0.z = Q[2][1]; c0.z = Q[2][2];

	// Think of it as #gipt different fiber orientations that are distributed within a plane of coordinate system (a0,b0)
	//  and that all contribute to the tangent in proportion with the probability of having a fiber in this orientation
	for(i=0;i<gipt;++i)
	{
	    // vector describing the fibers which make an angle pp[i] with vector a0 in plane (a0,b0)
		vec3d n0;
		n0.x = cos(pp[i])*a0.x+sin(pp[i])*b0.x;
		n0.y = cos(pp[i])*a0.y+sin(pp[i])*b0.y;
		n0.z = cos(pp[i])*a0.z+sin(pp[i])*b0.z;
		
		// probability of having a fiber along this vector 
		if (vmc==1) // Semi-circular von Mises distribution
		{
			distribution = 1. / ( pi * bessi0(kf) ) * exp( kf * cos( 2.* (pp[i]-tp) ) );
		}
		else if (vmc==2) // Constrained von Mises distribution
		{
			beta = pow(bessi1(kf)/bessi0(kf),var_n);  // Weighting factor (isotropic vs anisotropic)
			distribution = (1.-beta)/pi + beta / ( pi * bessi0(kf) ) * exp( kf * cos( 2.* (pp[i]-tp) ) );
		}		
		
	// add the fiber stress by Gauss integration : m_fib.Tangent(mp) is the contribution to the tangent of the fibers in direction pp[i], 
	// distribution is the probability of having a fiber in this direction, and wmg[i] is the Gauss Weight for integration
		c += wmg[i]*distribution*m_fib.DevFiberTangent(mp, n0);
	}

	return c;
}
