#pragma once
#include "FEElasticMaterial.h"

//-----------------------------------------------------------------------------
class FEJ2PlasticMaterialPoint : public FEMaterialPoint
{
public:
	FEJ2PlasticMaterialPoint(FEMaterialPoint* pt) : FEMaterialPoint(pt){}

	FEMaterialPoint* Copy()
	{
		FEJ2PlasticMaterialPoint* pt = new FEJ2PlasticMaterialPoint(*this);
		if (m_pNext) pt->m_pNext = m_pNext->Copy();
		return pt;
	}

	void Init()
	{
		// intialize data to zero
		e0.zero();
		e1.zero();
		sn.zero();
		b = false;
		Y1 = Y0;

		// don't forget to intialize the nested data
		FEMaterialPoint::Init();
	}

	void Update(const FETimeInfo& timeInfo)
	{
		FEElasticMaterialPoint& pt = *m_pNext->ExtractData<FEElasticMaterialPoint>();
		e0 = e1;
		sn = pt.m_s;
		Y0 = Y1;

		// don't forget to call the base class
		FEMaterialPoint::Update(timeInfo);
	}

	void Serialize(DumpStream& ar)
	{
		if (ar.IsSaving())
		{
			ar << e0 << e1 << sn;
			ar << Y0 << Y1 << b;
		}
		else
		{
			ar >> e0 >> e1 >> sn;
			ar >> Y0 >> Y1 >> b;
		}
		FEMaterialPoint::Serialize(ar);
	}

public:
	mat3ds	e0, e1;		// strain at time n and n+1
	mat3ds	sn;			// stress at time n
	double	Y0, Y1;		// yield strenght at time n, n+1
	bool	b;			// plasticity flag
};

//-----------------------------------------------------------------------------
//! This class implements a simple von-Mises plasticity model with isotropic
//! hardening. 
class FEVonMisesPlasticity : public FESolidMaterial
{
public:
	FEVonMisesPlasticity(FEModel* pfem);

public:
	double	m_E;	//!< Young's modulus
	double	m_v;	//!< Poisson's ratio

	double	m_K;	//!< bulk modulus
	double	m_G;	//!< shear modulus
	double	m_Y;	//!< initial yield strength
	double	m_H;	//!< hardening modulus 

public:
	FEMaterialPoint* CreateMaterialPointData() override;

	//! calculate stress at material point
	mat3ds Stress(FEMaterialPoint& pt) override;

	//! calculate tangent stiffness at material point
	tens4ds Tangent(FEMaterialPoint& pt) override;

	//! data initialization and checking
	bool Init() override;

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
