#pragma once
#include "FEMaterial.h"

//-----------------------------------------------------------------------------
//! This class implements a material that consists of a continuous ellipsoidal fiber distribution
//! superposed on a charged (swelling) gel described by the equations of Donnan equilibrium

//! This material is orignally due to Gerard Ateshian and is used to model
//! articular cartilage.

class FERandomFiberDonnanEquilibrium :	public FEElasticMaterial
	{
	public:
		FERandomFiberDonnanEquilibrium(void);
		
	public:
		//! calculate stress at material point
		virtual mat3ds Stress(FEMaterialPoint& pt);
		
		//! calculate tangent stiffness at material point
		virtual tens4ds Tangent(FEMaterialPoint& pt);
		
		//! material parameter intialization and checking
		void Init();

		//! return bulk modulus
		virtual double BulkModulus();
		
	public:
		double	m_phiwr;	// fluid volume fraction in reference configuration
		double	m_cFr;		// fixed charge density in reference configuration
		double	m_bosm;		// bath osmolarity
		double	m_Rgas;		// universal gas constant
		double	m_Tabs;		// absolute temperature
		
		double	m_beta[3];	// power in power-law relation
		double	m_ksi[3];	// coefficient in power-law relation
		
		// -------------------------------
		
		static	int	m_nres;	// integration rule
		
		static double	m_cth[];
		static double	m_sth[];
		static double	m_cph[];
		static double	m_sph[];
		static double	m_w[];
		
		// declare as registered
		DECLARE_REGISTERED(FERandomFiberDonnanEquilibrium);
		
		// declare the parameter list
		DECLARE_PARAMETER_LIST();
	};
