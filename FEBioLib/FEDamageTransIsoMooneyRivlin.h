#pragma once
#include "FEUncoupledMaterial.h"

//-----------------------------------------------------------------------------
// We first define a material point that stores the damage variable.
class FETIMRDamageMaterialPoint : public FEMaterialPoint
{
public:
	FETIMRDamageMaterialPoint(FEMaterialPoint *pt) : FEMaterialPoint(pt) {}

	FEMaterialPoint* Copy()
	{
		FETIMRDamageMaterialPoint* pt = new FETIMRDamageMaterialPoint(*this);
		if (m_pt) pt->m_pt = m_pt->Copy();
		return pt;
	}

	void Init(bool bflag)
	{
		FEElasticMaterialPoint& pt = *m_pt->ExtractData<FEElasticMaterialPoint>();
		if (bflag)
		{
			// intialize data to zero
			m_MEmax = 0;
			m_MEtrial = 0;
			m_Dm = 1;

			m_FEmax = 0;
			m_FEtrial = 0;
			m_Df = 1;
		}
		else
		{
			m_MEmax = max(m_MEmax, m_MEtrial);
			m_FEmax = max(m_FEmax, m_FEtrial);
		}

		// don't forget to intialize the nested data
		if (m_pt) m_pt->Init(bflag);
	}

	void Serialize(DumpFile& ar)
	{
		if (ar.IsSaving())
		{
			ar << m_MEmax << m_FEmax;
		}
		else
		{
			ar >> m_MEmax >> m_FEmax;
		}
	}

public:
	// matrix
	double	m_MEtrial;			//!< trial strain at time t
	double	m_MEmax;			//!< max strain variable up to time t
	double	m_Dm;				//!< damage reduction factor

	// fiber
	double	m_FEtrial;			//!< trial strain at time t
	double	m_FEmax;			//!< max strain variable up to time t
	double	m_Df;				//!< damage reduction factor
};

//-----------------------------------------------------------------------------
class FEDamageTransIsoMooneyRivlin : public FEUncoupledMaterial
{
public:
	FEDamageTransIsoMooneyRivlin(void);

public:
	// Mooney-Rivlin parameters
	double	m_c1;	//!< Mooney-Rivlin coefficient C1
	double	m_c2;	//!< Mooney-Rivlin coefficient C2

	// fiber parameters
	double	m_c3;
	double	m_c4;

	// Matrix damage parameters
	double	m_Mbeta;		//!< damage parameter beta
	double	m_Msmin;		//!< damage parameter psi-min
	double	m_Msmax;		//!< damage parameter psi-max

	// Fiber damage parameters
	double	m_Fbeta;
	double	m_Fsmin;
	double	m_Fsmax;

public:
	// returns a pointer to a new material point object
	virtual FEMaterialPoint* CreateMaterialPointData() { return new FETIMRDamageMaterialPoint(new FEElasticMaterialPoint); }

public:
	//! calculate deviatoric stress at material point
	mat3ds DevStress(FEMaterialPoint& pt);

	//! calculate deviatoric tangent stiffness at material point
	tens4ds DevTangent(FEMaterialPoint& pt);

	//! data initialization
	void Init();

protected:
	mat3ds MatrixStress(FEMaterialPoint& mp);
	mat3ds FiberStress (FEMaterialPoint& mp);
	tens4ds MatrixTangent(FEMaterialPoint& pt);
	tens4ds FiberTangent (FEMaterialPoint& pt);

protected:
	// calculate damage reduction factor for matrix
	double MatrixDamage(FEMaterialPoint& pt);

	// calculate damage reduction factor for fibers
	double FiberDamage(FEMaterialPoint& pt);

	double MatrixDamageDerive(FEMaterialPoint& pt);
	double FiberDamageDerive(FEMaterialPoint& pt);

public:

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
