/*
 *  FEEllipsoidalFiberDistribution.h
 *  FEBioXCode
 *
 *  Created by Gerard Ateshian on 2/17/10.
 *  Copyright 2010 Columbia University. All rights reserved.
 *
 */

#include "FEMaterial.h"

void FEEllipsoidalFiberDistributionInit(const double m_ksi[], const double m_beta[]);

mat3ds FEEllipsoidalFiberDistributionStress(const double m_ksi[], const double m_beta[], 
											const int nint, const double m_cth[], const double m_sth[], 
											const double m_cph[], const double m_sph[], const double m_w[], 
											const double J, const mat3d F, const mat3d Q);

tens4ds FEEllipsoidalFiberDistributionTangent(const double m_ksi[], const double m_beta[], 
											  const int nint, const double m_cth[], const double m_sth[], 
											  const double m_cph[], const double m_sph[], const double m_w[], 
											  const double J, const mat3d F, const mat3d Q);

