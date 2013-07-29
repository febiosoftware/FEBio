#pragma once
#include "FEElasticMaterial.h"

//-----------------------------------------------------------------------------
// This material is a first attempt to include damage in hyper-elastic materials.
// It assumes the simple damage model as defined in Simo, CMAME 60 (1987), 153-173

//-----------------------------------------------------------------------------
// We first define a material point that stores the damage variable.
class FEDamageMaterialPoint : public FEMaterialPoint
{
public:
	FEDamageMaterialPoint(FEMaterialPoint *pt) : FEMaterialPoint(pt) {}

	FEMaterialPoint* Copy()
	{
		FEDamageMaterialPoint* pt = new FEDamageMaterialPoint(*this);
		if (m_pt) pt->m_pt = m_pt->Copy();
		return pt;
	}

	void Init(bool bflag)
	{
		FEElasticMaterialPoint& pt = *m_pt->ExtractData<FEElasticMaterialPoint>();
		if (bflag)
		{
			// intialize data to zero
			m_Emax = 0;
			m_Etrial = 0;
			m_D = 1;
		}
		else
		{
			m_Emax = max(m_Emax, m_Etrial);
		}

		// don't forget to intialize the nested data
		if (m_pt) m_pt->Init(bflag);
	}

	void Serialize(DumpFile& ar)
	{
		if (ar.IsSaving())
		{
			ar << m_Emax;
		}
		else
		{
			ar >> m_Emax;
		}
	}

public:
	double	m_Etrial;		//!< trial strain at time t
	double	m_Emax;			//!< max strain variable up to time t
	double	m_D;			//!< damage reduction factor
};

//-----------------------------------------------------------------------------
class FEDamageNeoHookean : public FEElasticMaterial
{
public:
	FEDamageNeoHookean(void);

public:
	double	m_E;	//!< Young's modulus
	double	m_v;	//!< Poisson's ratio

	double	m_alpha;	//!< damage parameter alpha
	double	m_beta;		//!< damage parameter beta

protected:
	double	m_lam;
	double	m_mu;

public:
	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt);

	//! calculate tangent stiffness at material point
	virtual tens4ds Tangent(FEMaterialPoint& pt);

	//! data initialization and checking
	void Init();

	// returns a pointer to a new material point object
	virtual FEMaterialPoint* CreateMaterialPointData() 
	{ 
		return new FEDamageMaterialPoint(new FEElasticMaterialPoint);
	}

	// calculate damage reduction factor
	double Damage(FEMaterialPoint& pt);

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
