/*
 *  FEDonnanEquilibrium.h
 *  FEBioXCode
 *
 *  Created by Gerard Ateshian on 2/17/10.
 *
 */

#include "FEMaterial.h"

//-----------------------------------------------------------------------------
//! Material class that implements Donnan equilibrium. 
//
class FEDonnanEquilibrium : public FEMaterial
{
public:
	//! Initialization routine
	void Init();

	//! Returns the Cauchy stress
	mat3ds Stress(FEMaterialPoint& mp);

	//! Returs the spatial tangent
	tens4ds Tangent(FEMaterialPoint& mp);

	//! returns the bulkmodulus
	double BulkModulus();

public:
	double	m_phiwr;	//!< fluid volume fraction in reference configuration
	double	m_cFr;		//!< fixed charge density in reference configuration
	double	m_Rgas;		//!< universal gas constant
	double	m_Tabs;		//!< absolute temperature
	double	m_bosm;		//!< bath osmolarity
};
