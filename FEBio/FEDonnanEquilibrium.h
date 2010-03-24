/*
 *  FEDonnanEquilibrium.h
 *  FEBioXCode
 *
 *  Created by Gerard Ateshian on 2/17/10.
 *
 */

#include "FEMaterial.h"

void FEDonnanEquilibriumInit(const double m_phiwr, const double m_cFr, 
							 const double m_Rgas, const double m_Tabs, 
							 const double m_bosm);

mat3ds FEDonnanEquilibriumStress(const double m_phiwr, const double m_cFr, 
								 const double m_Rgas, const double m_Tabs, 
								 const double m_bosm, const double J);

tens4ds FEDonnanEquilibriumTangent(const double m_phiwr, const double m_cFr, 
								   const double m_Rgas, const double m_Tabs, 
								   const double m_bosm, const double J);

double FEDonnanEquilibriumBulkModulus(const double m_phiwr, const double m_cFr, 
							 const double m_Rgas, const double m_Tabs, 
							 const double m_bosm);

