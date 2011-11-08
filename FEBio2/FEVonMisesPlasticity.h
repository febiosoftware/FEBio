#pragma once
#include "FECore/FEElasticMaterial.h"

//-----------------------------------------------------------------------------
class FEJ2PlasticMaterialPoint : public FEMaterialPoint
{
public:
	FEJ2PlasticMaterialPoint(FEMaterialPoint* pt) : FEMaterialPoint(pt){}

	FEMaterialPoint* Copy()
	{
		FEJ2PlasticMaterialPoint* pt = new FEJ2PlasticMaterialPoint(*this);
		if (m_pt) pt->m_pt = m_pt->Copy();
		return pt;
	}

	void Init(bool bflag)
	{
		FEElasticMaterialPoint& pt = *m_pt->ExtractData<FEElasticMaterialPoint>();
		if (bflag)
		{
			// intialize data to zero
			e0.zero();
			e1.zero();
			sn.zero();
			b = false;
		}
		else
		{
			e0 = e1;
			sn = pt.s;
		}

		// don't forget to intialize the nested data
		if (m_pt) m_pt->Init(bflag);
	}

	void Serialize(DumpFile& ar)
	{
		if (ar.IsSaving())
		{
		}
		else
		{
		}
	}

public:
	mat3ds	e0, e1;	// strain at time n and n+1
	mat3ds	sn;		// stress at time n
	bool	b;		// plasticity flag
};

//-----------------------------------------------------------------------------
class FEVonMisesPlasticity : public FEElasticMaterial
{
public:
	FEVonMisesPlasticity(void);

public:
	double	m_K;	// bulk modulus
	double	m_G;	// shear modulus
	double	m_Y;	// Yield-stress

public:
	virtual FEMaterialPoint* CreateMaterialPointData() { return new FEJ2PlasticMaterialPoint(new FEElasticMaterialPoint); }

	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt);

	//! calculate tangent stiffness at material point
	virtual tens4ds Tangent(FEMaterialPoint& pt);

	//! data initialization and checking
	void Init();

	// declare as registered
	DECLARE_REGISTERED(FEVonMisesPlasticity);

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
