#pragma once
#include "FEUncoupledMaterial.h"

//-----------------------------------------------------------------------------
//! Uncoupled elastic mixtures

//! This class describes a mixture of uncoupled elastic solids.  The user must declare
//! uncoupled elastic solids that can be combined within this class.  The stress and
//! tangent tensors evaluated in this class represent the sum of the respective
//! tensors of all the solids forming the mixture.

class FEUncoupledElasticMixture : public FEUncoupledMaterial
	{
	public:
		FEUncoupledElasticMixture() {}
		
	public:
		vector <FEUncoupledMaterial*>	m_pMat;	//!< pointers to elastic materials
		
	public:
		//! calculate stress at material point
		mat3ds DevStress(FEMaterialPoint& pt);
		
		//! calculate tangent stiffness at material point
		tens4ds DevTangent(FEMaterialPoint& pt);
		
		//! data initialization and checking
		void Init();
		
		// declare as registered
		DECLARE_REGISTERED(FEUncoupledElasticMixture);
		
		// declare the parameter list
//		DECLARE_PARAMETER_LIST();
	};
