#pragma once
#include "FEBiphasic.h"

//-----------------------------------------------------------------------------
// This class implements a poroelastic material that has a strain-dependent
// permeability which is orthotropic in the reference state, but exhibits
// further strain-induced anisotropy, according to the constitutive relation
// of Ateshian and Weiss (JBME 2010)

class FEPermRefTransIso :	public FEHydraulicPermeability
	{
	public:
		//! constructor
		FEPermRefTransIso(FEModel* pfem);
		
		//! permeability
		mat3ds Permeability(FEMaterialPoint& pt);
		
		//! Tangent of permeability
		tens4ds Tangent_Permeability_Strain(FEMaterialPoint& mp);
		
	public:
		double	m_perm0;		//!< permeability for I term
		double	m_perm1T;		//!< transverse permeability for b term
		double	m_perm1A;		//!< axial permeability for b term
		double	m_perm2T;		//!< transverse permeability for b^2 term
		double	m_perm2A;		//!< axial permeability for b^2 term
		double	m_M0;			//!< nonlinear exponential coefficient for I term
		double	m_MT;			//!< nonlinear exponential coefficient for transverse direction
		double	m_MA;			//!< nonlinear exponential coefficient for axial direction
		double	m_alpha0;		//!< nonlinear power exponent for I term
		double	m_alphaT;		//!< nonlinear power exponent for transverse direction
		double	m_alphaA;		//!< nonlinear power exponent for axial direction
		
		// declare parameter list
		DECLARE_PARAMETER_LIST();
	};
