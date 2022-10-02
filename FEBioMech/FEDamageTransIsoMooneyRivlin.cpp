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
#include "FEDamageTransIsoMooneyRivlin.h"

FETIMRDamageMaterialPoint::FETIMRDamageMaterialPoint(FEMaterialPointData*pt) : FEMaterialPointData(pt) {}

FEMaterialPointData* FETIMRDamageMaterialPoint::Copy()
{
	FETIMRDamageMaterialPoint* pt = new FETIMRDamageMaterialPoint(*this);
	if (m_pNext) pt->m_pNext = m_pNext->Copy();
	return pt;
}

void FETIMRDamageMaterialPoint::Init()
{
	FEMaterialPointData::Init();

	// intialize data to zero
	m_MEmax = 0;
	m_MEtrial = 0;
	m_Dm = 0;

	m_FEmax = 0;
	m_FEtrial = 0;
	m_Df = 0;
}

void FETIMRDamageMaterialPoint::Update(const FETimeInfo& timeInfo)
{
	FEMaterialPointData::Update(timeInfo);

	m_MEmax = max(m_MEmax, m_MEtrial);
	m_FEmax = max(m_FEmax, m_FEtrial);
}

void FETIMRDamageMaterialPoint::Serialize(DumpStream& ar)
{
	FEMaterialPointData::Serialize(ar);
	ar & m_MEtrial & m_MEmax & m_Dm;
	ar & m_FEtrial & m_FEmax & m_Df;
}


// define the material parameters
BEGIN_FECORE_CLASS(FEDamageTransIsoMooneyRivlin, FEUncoupledMaterial)
	ADD_PARAMETER(m_c1, FE_RANGE_GREATER(0.0), "c1");
	ADD_PARAMETER(m_c2, "c2");
	ADD_PARAMETER(m_c3, "c3");
	ADD_PARAMETER(m_c4, FE_RANGE_GREATER(0.0), "c4");
	ADD_PARAMETER(m_Mbeta, "Mbeta");
	ADD_PARAMETER(m_Msmin, "Msmin");
	ADD_PARAMETER(m_Msmax, "Msmax");
	ADD_PARAMETER(m_Fbeta, "Fbeta");
	ADD_PARAMETER(m_Fsmin, "Fsmin");
	ADD_PARAMETER(m_Fsmax, "Fsmax");

	ADD_PROPERTY(m_Q, "mat_axis")->SetFlags(FEProperty::Optional);

END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEDamageTransIsoMooneyRivlin::FEDamageTransIsoMooneyRivlin(FEModel* pfem) : FEUncoupledMaterial(pfem)
{

}

//-----------------------------------------------------------------------------
FEMaterialPointData* FEDamageTransIsoMooneyRivlin::CreateMaterialPointData() 
{ 
	return new FETIMRDamageMaterialPoint(new FEElasticMaterialPoint); 
}

//-----------------------------------------------------------------------------
mat3ds FEDamageTransIsoMooneyRivlin::DevStress(FEMaterialPoint& mp)
{
	// matrix stress
	mat3ds sm = MatrixStress(mp);

	// matrix reduction factor
	double gm = MatrixDamage(mp);

	// fiber stress
	mat3ds sf = FiberStress(mp);

	// fiber reduction factor
	double gf = FiberDamage(mp);

	return sm*gm + sf*gf;
}

//-----------------------------------------------------------------------------
//! Calculate the deviatoric stress
mat3ds FEDamageTransIsoMooneyRivlin::MatrixStress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// determinant of deformation gradient
	double J = pt.m_J;

	// calculate deviatoric left Cauchy-Green tensor
	mat3ds B = pt.DevLeftCauchyGreen();

	// calculate square of B
	mat3ds B2 = B.sqr();

	// Invariants of B (= invariants of C)
	// Note that these are the invariants of Btilde, not of B!
	double I1 = B.tr();

	// --- TODO: put strain energy derivatives here ---
	//
	// W = C1*(I1 - 3) + C2*(I2 - 3)
	//
	// Wi = dW/dIi
	double W1 = m_c1;
	double W2 = m_c2;
	// ---

	// calculate T = F*dW/dC*Ft
	// T = F*dW/dC*Ft
	mat3ds T = B*(W1 + W2*I1) - B2*W2;

	return T.dev()*(2.0/J);
}

//-----------------------------------------------------------------------------
// Calculate the fiber stress
mat3ds FEDamageTransIsoMooneyRivlin::FiberStress(FEMaterialPoint &mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// deformation gradient
	mat3d& F = pt.m_F;
	double J = pt.m_J;
	double Jm13 = pow(J, -1.0/3.0);

	// get the local coordinate systems
	mat3d Q = GetLocalCS(mp);

	// get the initial fiber direction
	vec3d a0 = Q.col(0);

	// calculate the current material axis lam*a = F*a0;
	vec3d a = F*a0;

	// normalize material axis and store fiber stretch
	double lam, lamd;
	lam = a.unit();
	lamd = lam*Jm13; // i.e. lambda tilde

	// invariant I4
	double I4 = lamd*lamd;

	// strain energy derivative
	double W4 = (I4 - 1)*m_c3*exp(m_c4*(I4-1)*(I4-1));

	// calculate dyad of a: AxA = (a x a)
	mat3ds AxA = dyad(a);

	// calculate FdWf/dCFt = I4*W4*(a x a)
	mat3ds T = AxA*(W4*I4);

	// return stress
	return T.dev()*(2.0/J);
}

//-----------------------------------------------------------------------------
tens4ds FEDamageTransIsoMooneyRivlin::DevTangent(FEMaterialPoint& mp)
{
	// matrix tangent
	tens4ds Cm = MatrixTangent(mp);

	// matrix damage
	double gm = MatrixDamage(mp);

	// fiber tangent
	tens4ds Cf = FiberTangent(mp);

	// fiber damage
	double gf = FiberDamage(mp);

	return Cm*gm + Cf*gf;
}

//-----------------------------------------------------------------------------
//! Calculate the deviatoric tangent
tens4ds FEDamageTransIsoMooneyRivlin::MatrixTangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// determinant of deformation gradient
	double J = pt.m_J;
	double Ji = 1.0/J;

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
	W1 = m_c1;
	W2 = m_c2;
	// ---

	// calculate dWdC:C
	double WC = W1*I1 + 2*W2*I2;

	// calculate C:d2WdCdC:C
	double CWWC = 2*I2*W2;

	// deviatoric cauchy-stress, trs = trace[s]/3
	mat3ds devs = pt.m_s.dev();

	// Identity tensor
	mat3ds I(1,1,1,0,0,0);

	tens4ds IxI = dyad1s(I);
	tens4ds I4  = dyad4s(I);
	tens4ds BxB = dyad1s(B);
	tens4ds B4  = dyad4s(B);

	// d2W/dCdC:C
	mat3ds WCCxC = B*(W2*I1) - B2*W2;

	tens4ds cw = (BxB - B4)*(W2*4.0*Ji) - dyad1s(WCCxC, I)*(4.0/3.0*Ji) + IxI*(4.0/9.0*Ji*CWWC);

	tens4ds c = dyad1s(devs, I)*(-2.0/3.0) + (I4 - IxI/3.0)*(4.0/3.0*Ji*WC) + cw;

	return c;
}


//-----------------------------------------------------------------------------
// Fiber material tangent
//
tens4ds FEDamageTransIsoMooneyRivlin::FiberTangent(FEMaterialPoint &mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// deformation gradient
	mat3d& F = pt.m_F;
	double J = pt.m_J;
	double Jm13 = pow(J, -1.0/3.0);
	double Ji = 1.0/J;

	// get the local coordinate systems
	mat3d Q = GetLocalCS(mp);

	// get initial local material axis
	vec3d a0 = Q.col(0);

	// calculate current local material axis
	vec3d a = F*a0;

	double lam = a.unit();

	// deviatoric stretch
	double lamd = lam*Jm13;

	double I4 = lamd*lamd;

	// strain energy derivative
	double W4 = (I4 - 1)*m_c3*exp(m_c4*(I4-1)*(I4-1));
	double W44 = m_c3*exp(m_c4*(I4 - 1)*(I4 - 1)) + 2*m_c3*m_c4*(I4 - 1)*(I4 - 1)*exp(m_c4*(I4 - 1)*(I4 - 1));

	// calculate dWdC:C
	double WC = W4*I4;

	// calculate C:d2WdCdC:C
	double CWWC = W44*I4*I4;

	mat3dd I(1);	// Identity
	tens4ds IxI = dyad1s(I);
	tens4ds Id4  = dyad4s(I);

	mat3ds AxA = dyad(a);
	tens4ds AxAxAxA = dyad1s(AxA);

	tens4ds cw = AxAxAxA*(4.0*Ji*W44*I4*I4) - dyad1s(I, AxA)*(4.0/3.0*Ji*W44*I4*I4);

	tens4ds c = (Id4 - IxI/3.0)*(4.0/3.0*Ji*WC) + IxI*(4.0/9.0*Ji*CWWC) + cw;

	// see if we need to add the stress
	FETIMRDamageMaterialPoint& dp = *mp.ExtractData<FETIMRDamageMaterialPoint>();
	if (dp.m_FEtrial > dp.m_FEmax)
	{
		mat3ds devs = pt.m_s.dev();
		double dg = FiberDamageDerive(mp);
		c += dyad1s(devs)*(J*dg/dp.m_FEtrial);
	}

	return c;
}

//-----------------------------------------------------------------------------
double FEDamageTransIsoMooneyRivlin::DevStrainEnergyDensity(FEMaterialPoint& mp)
{
	// matrix sed
	double sedm = MatrixStrainEnergyDensity(mp);
    
	// matrix reduction factor
	double gm = MatrixDamage(mp);
    
	// fiber sed
	double sedf = FiberStrainEnergyDensity(mp);
    
	// fiber reduction factor
	double gf = FiberDamage(mp);
    
	return sedm*gm + sedf*gf;
}

//-----------------------------------------------------------------------------
//! Calculate the matrix strain energy density
double FEDamageTransIsoMooneyRivlin::MatrixStrainEnergyDensity(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
	// calculate deviatoric left Cauchy-Green tensor
	mat3ds B = pt.DevLeftCauchyGreen();
    
	// calculate square of B
	mat3ds B2 = B.sqr();
    
	// Invariants of B (= invariants of C)
	// Note that these are the invariants of Btilde, not of B!
	double I1 = B.tr();
	double I2 = 0.5*(I1*I1 - B2.tr());
    
	// --- TODO: put strain energy derivatives here ---
	//
	// W = C1*(I1 - 3) + C2*(I2 - 3)
	//
    double sed = m_c1*(I1 - 3) + m_c2*(I2 - 3);
    
	return sed;
}

//-----------------------------------------------------------------------------
// Calculate the fiber strain energy density
double FEDamageTransIsoMooneyRivlin::FiberStrainEnergyDensity(FEMaterialPoint &mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
    
	// deformation gradient
	mat3d& F = pt.m_F;
	double J = pt.m_J;
	double Jm13 = pow(J, -1.0/3.0);

	// get the local coordinate systems
	mat3d Q = GetLocalCS(mp);

	// get the initial fiber direction
	vec3d a0 = Q.col(0);
    
	// calculate the current material axis lam*a = F*a0;
	vec3d a = F*a0;
    
	// normalize material axis and store fiber stretch
	double lam, lamd;
	lam = a.unit();
	lamd = lam*Jm13; // i.e. lambda tilde
    
	// invariant I4
	double I4 = lamd*lamd;
    
	// strain energy derivative
	double sed = 0.5*m_c3/m_c4*(exp(m_c4*(I4-1)*(I4-1))-1);
    
	// return sed
	return sed;
}

//-----------------------------------------------------------------------------
// Calculate damage reduction factor for matrix
double FEDamageTransIsoMooneyRivlin::MatrixDamage(FEMaterialPoint &mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// calculate right Cauchy-Green tensor
	mat3ds C = pt.DevRightCauchyGreen();
	mat3ds C2 = C.sqr();

	// Invariants
	double I1 = C.tr();
	double I2 = 0.5*(I1*I1 - C2.tr());

	// strain-energy value
	double SEF = m_c1*(I1 - 3) + m_c2*(I2 - 3);

	// get the damage material point data
	FETIMRDamageMaterialPoint& dp = *mp.ExtractData<FETIMRDamageMaterialPoint>();

	// calculate trial-damage parameter
	dp.m_MEtrial = sqrt(2.0*fabs(SEF));

	// calculate damage parameter
	double Es = max(dp.m_MEtrial, dp.m_MEmax);

	// calculate reduction parameter
	double g = 1.0;
	if (Es < m_Msmin) g = 1.0;
	else if (Es > m_Msmax) g = 0.0;
	else 
	{
		double F = (Es - m_Msmin)/(m_Msmin - m_Msmax);
		g = 1.0 - (1.0 - m_Mbeta + m_Mbeta*F*F)*(F*F);
	}

	dp.m_Dm = 1-g;
	return g;
}

//-----------------------------------------------------------------------------
// Calculate damage reduction factor for matrix
double FEDamageTransIsoMooneyRivlin::MatrixDamageDerive(FEMaterialPoint &mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// calculate right Cauchy-Green tensor
	mat3ds C = pt.DevRightCauchyGreen();
	mat3ds C2 = C.sqr();

	// Invariants
	double I1 = C.tr();
	double I2 = 0.5*(I1*I1 - C2.tr());

	// strain-energy value
	double SEF = m_c1*(I1 - 3) + m_c2*(I2 - 3);

	// get the damage material point data
	FETIMRDamageMaterialPoint& dp = *mp.ExtractData<FETIMRDamageMaterialPoint>();

	// calculate trial-damage parameter
	dp.m_MEtrial = sqrt(2.0*fabs(SEF));

	// calculate damage parameter
	double Es = max(dp.m_MEtrial, dp.m_MEmax);

	// calculate reduction parameter
	double dg = 0.0;
	if (Es < m_Msmin) dg = 0.0;
	else if (Es > m_Msmax) dg = 0.0;
	else 
	{
		double h = m_Msmax - m_Msmin;
		double F = (Es - m_Msmin)/h;
		dg = -2.0*F/h*(1 - m_Mbeta + m_Mbeta*F*F)-(F*F)*(2.0*m_Mbeta*F/h);
	}

	return dg;
}


//-----------------------------------------------------------------------------
// Calculate damage reduction factor for fibers
double FEDamageTransIsoMooneyRivlin::FiberDamage(FEMaterialPoint &mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// deformation gradient
	mat3d& F = pt.m_F;
	double J = pt.m_J;
	double Jm13 = pow(J, -1.0/3.0);

	// get the local coordinate systems
	mat3d Q = GetLocalCS(mp);

	// get the initial fiber direction
	vec3d a0 = Q.col(0);

	// calculate the current material axis lam*a = F*a0;
	vec3d a = F*a0;

	// normalize material axis and store fiber stretch
	double lam, lamd;
	lam = a.unit();
	lamd = lam*Jm13; // i.e. lambda tilde

	// invariant I4
	double I4 = lamd*lamd;

	// strain energy value
	double SEF = 0.5*m_c3/m_c4*(exp(m_c4*(I4-1)*(I4-1))-1);

	// get the damage material point data
	FETIMRDamageMaterialPoint& dp = *mp.ExtractData<FETIMRDamageMaterialPoint>();

	// calculate trial-damage parameter
	dp.m_FEtrial = sqrt(2.0*fabs(SEF));

	// calculate damage parameter
	double Es = max(dp.m_FEtrial, dp.m_FEmax);

	// calculate reduction parameter
	double g = 1.0;
	if (Es < m_Fsmin) g = 1.0;
	else if (Es > m_Fsmax) g = 0.0;
	else 
	{
		double F = (Es - m_Fsmin)/(m_Fsmin - m_Fsmax);
		g = 1.0 - (1.0 - m_Fbeta + m_Fbeta*F*F)*(F*F);
	}

	dp.m_Df = 1-g;
	return g;
}


//-----------------------------------------------------------------------------
// Calculate damage reduction factor for fibers
double FEDamageTransIsoMooneyRivlin::FiberDamageDerive(FEMaterialPoint &mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// deformation gradient
	mat3d& F = pt.m_F;
	double J = pt.m_J;
	double Jm13 = pow(J, -1.0/3.0);

	// get the local coordinate systems
	mat3d Q = GetLocalCS(mp);

	// get the initial fiber direction
	vec3d a0 = Q.col(0);

	// calculate the current material axis lam*a = F*a0;
	vec3d a = F*a0;

	// normalize material axis and store fiber stretch
	double lam, lamd;
	lam = a.unit();
	lamd = lam*Jm13; // i.e. lambda tilde

	// invariant I4
	double I4 = lamd*lamd;

	// strain energy value
	double SEF = 0.5*m_c3/m_c4*(exp(m_c4*(I4-1)*(I4-1))-1);

	// get the damage material point data
	FETIMRDamageMaterialPoint& dp = *mp.ExtractData<FETIMRDamageMaterialPoint>();

	// calculate trial-damage parameter
	dp.m_FEtrial = sqrt(2.0*fabs(SEF));

	// calculate damage parameter
	double Es = max(dp.m_FEtrial, dp.m_FEmax);

	// calculate reduction parameter
	double dg = 0.0;
	if (Es < m_Fsmin) dg = 0.0;
	else if (Es > m_Fsmax) dg = 0.0;
	else 
	{
		double h = m_Fsmin - m_Fsmax;
		double F = (Es - m_Fsmin)/h;
		dg = -2.0*F/h*(1 - m_Fbeta + m_Fbeta*F*F)-(F*F)*(2.0*m_Fbeta*F/h);
	}

	return dg;
}

//-----------------------------------------------------------------------------
//! damage
double FEDamageTransIsoMooneyRivlin::Damage(FEMaterialPoint& pt)
{
    return MatrixDamage(pt) + FiberDamage(pt);
}
