#include "FEUncoupledMaterial.h"

//-----------------------------------------------------------------------------
//! Material class for the uncoupled ellipsoidal fiber distribution
//!

class FEEFDUncoupled : public FEUncoupledMaterial
	{
	public:
		FEEFDUncoupled() {m_unstable = true;}
		
		//! Initialization
		void Init();
		
		//! Cauchy stress
		virtual mat3ds DevStress(FEMaterialPoint& mp);
		
		// Spatial tangent
		virtual tens4ds DevTangent(FEMaterialPoint& mp);
		
		// declare as registered
		DECLARE_REGISTERED(FEEFDUncoupled);
		
		// declare the parameter list
		DECLARE_PARAMETER_LIST();
		
	public:
		double	m_beta[3];	// power in power-law relation
		double	m_ksi[3];	// coefficient in power-law relation
		
		static int		m_nres;	// integration rule
		static double	m_cth[];
		static double	m_sth[];
		static double	m_cph[];
		static double	m_sph[];
		static double	m_w[];
	};
