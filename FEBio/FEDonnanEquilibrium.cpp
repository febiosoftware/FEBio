/*
 *  FEDonnanEquilibrium.cpp
 *  FEBioXCode
 *
 *  Created by Gerard Ateshian on 2/17/10.
 *
 */

#include "FEDonnanEquilibrium.h"

//////////////////////////////////////////////////////////////////////
// FEDonnanEquilibrium
//////////////////////////////////////////////////////////////////////

void FEDonnanEquilibriumInit(const double m_phiwr, const double m_cFr, 
							 const double m_Rgas, const double m_Tabs, 
							 const double m_bosm)
{
	if (m_phiwr < 0 || m_phiwr > 1) throw MaterialError("phiw0 must be between 0. and 1.");
	if (m_Rgas < 0) throw MaterialError("R must be positive.");
	if (m_Tabs < 0) throw MaterialError("T must be positive.");
	if (m_bosm < 0) throw MaterialError("bosm must be positive.");
}

mat3ds FEDonnanEquilibriumStress(const double m_phiwr, const double m_cFr, 
								 const double m_Rgas, const double m_Tabs, 
								 const double m_bosm, const double J)
{
	
	// calculate fixed charge density in current configuration
	double cF = m_phiwr*m_cFr/(J-1+m_phiwr);
	
	// calculate osmotic pressure
	double p = m_Rgas*m_Tabs*(sqrt(cF*cF+m_bosm*m_bosm) - m_bosm);
	
	// calculate T = -p*I
	mat3dd I(1.0);	// identity tensor
	mat3ds s = -p*I;
	return s;
}

tens4ds FEDonnanEquilibriumTangent(const double m_phiwr, const double m_cFr, 
								   const double m_Rgas, const double m_Tabs, 
								   const double m_bosm, const double J)
{

	// calculate fixed charge density in current configuration
	double cF = m_phiwr*m_cFr/(J-1+m_phiwr);
	
	// calculate osmotic pressure
	double tosm = sqrt(cF*cF+m_bosm*m_bosm);	// tissue osmolarity
	double p = m_Rgas*m_Tabs*(tosm - m_bosm);	// osmotic pressure
	
	// calculate derivative of osmotic pressure w.r.t. J
	double bpi = m_Rgas*m_Tabs*J*cF*cF/(J-1+m_phiwr)/tosm;
	
	mat3dd I(1.0);	// Identity
	
	tens4ds IxI = dyad1s(I);
	tens4ds I4  = dyad4s(I);
	
	// calculate tangent osmotic modulus
	tens4ds c = bpi*IxI + p*(2.0*I4 - IxI);
	return c;
}

double FEDonnanEquilibriumBulkModulus(const double m_phiwr, const double m_cFr, 
									  const double m_Rgas, const double m_Tabs, 
									  const double m_bosm)
{
	// calculate osmotic pressure (assume J=1)
	double tosm = sqrt(m_cFr*m_cFr+m_bosm*m_bosm);	// tissue osmolarity
	double pi = m_Rgas*m_Tabs*(tosm - m_bosm);	// osmotic pressure
	
	// calculate derivative of osmotic pressure w.r.t. J
	double bpi = m_Rgas*m_Tabs*m_cFr*m_cFr/m_phiwr/tosm;
	
	return -(pi/3.0 + bpi);
}

