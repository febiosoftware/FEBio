#pragma once
#include "FECore/FEMaterial.h"

//-----------------------------------------------------------------------------
//! Material point data for visco-elastic materials
class FEViscoElasticMaterialPoint : public FEMaterialPoint
{
public:
	enum { MAX_TERMS = 6 };

public:
	FEViscoElasticMaterialPoint(FEMaterialPoint *pt) : FEMaterialPoint(pt) {}

	FEMaterialPoint* Copy()
	{
		FEViscoElasticMaterialPoint* pt = new FEViscoElasticMaterialPoint(*this);
		if (m_pt) pt->m_pt = m_pt->Copy();
		return pt;
	}

	void Init(bool bflag)
	{
		FEElasticMaterialPoint& pt = *m_pt->ExtractData<FEElasticMaterialPoint>();
		if (bflag)
		{
			// intialize data to zero
			m_se.zero();
			m_Sep.zero();
			for (int i=0; i<MAX_TERMS; ++i) { m_H[i].zero(); m_Hp[i].zero(); };
		}
		else
		{
			// the elastic stress stored in pt is the Cauchy stress.
			// however, we need to store the 2nd PK stress
			m_Sep = pt.pull_back(m_se);

			// copy previous data
			for (int i=0; i<MAX_TERMS; ++i) m_Hp[i] = m_H[i];
		}

		// don't forget to intialize the nested data
		if (m_pt) m_pt->Init(bflag);
	}

	void Serialize(DumpFile& ar)
	{
		if (m_pt) m_pt->Serialize(ar);

		if (ar.IsSaving())
		{
			ar << m_se;
			ar << m_Sep;
			ar << (int) MAX_TERMS;
			for (int i=0; i<MAX_TERMS; ++i) ar << m_H[i] << m_Hp[i];
		}
		else
		{
			ar >> m_se;
			ar >> m_Sep;
			int n;
			ar >> n;
			assert(n == MAX_TERMS);
			for (int i=0; i<MAX_TERMS; ++i) ar >> m_H[i] >> m_Hp[i];
		}
	}

public:
	mat3ds	m_se;	//!< elastic Cauchy stress
	mat3ds	m_Sep;	//!< elastic 2nd PK stress at previous time

	mat3ds	m_H[MAX_TERMS];		//!< internal variables
	mat3ds	m_Hp[MAX_TERMS];	//!< internal variables at previous timestep
};


//-----------------------------------------------------------------------------
//! This class implements a large deformation visco-elastic material

class FEViscoElasticMaterial :	public FENestedMaterial
{
public:
	// NOTE: make sure that this parameter is the 
	//       same as the MAX_TERMS in the FEViscoElasticMaterialPoint class
	enum { MAX_TERMS = FEViscoElasticMaterialPoint::MAX_TERMS };

public:
	//! default constructor
	FEViscoElasticMaterial();

	//! stress function
	mat3ds Stress(FEMaterialPoint& pt);

	//! tangent function
	virtual tens4ds Tangent(FEMaterialPoint& pt);

	// returns a pointer to a new material point object
	virtual FEMaterialPoint* CreateMaterialPointData() 
	{ 
		return new FEViscoElasticMaterialPoint(m_pBase->CreateMaterialPointData());
	}

public:
	double	m_t[MAX_TERMS];	//!< relaxation times
	double	m_g0;			//!< intitial visco-elastic coefficient
	double	m_g[MAX_TERMS];	//!< visco-elastic coefficients

public:
	// declare as registered
	DECLARE_REGISTERED(FEViscoElasticMaterial);

	// declare parameter list
	DECLARE_PARAMETER_LIST();
};
