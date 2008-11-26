#pragma once

#include "mat3d.h"
#include "Archive.h"

//-----------------------------------------------------------------------------
//! Material point class

//! This class implements the concept of a material point. This point carries
//! with it not only information about its location, both in the reference and  
//! current configuration but also about the local deformation. In addition
//! it contains the state information that is associated with the current
//! point.

class FEMaterialPoint
{
public:
	FEMaterialPoint(FEMaterialPoint* ppt = 0) : m_pt(ppt){}
	virtual ~FEMaterialPoint() { if (m_pt) { delete m_pt; m_pt = 0; } }

	//! The init function is used to intialize data
	virtual void Init(bool bflag) = 0;

	virtual FEMaterialPoint* Copy() = 0;

	virtual void Serialize(Archive& ar) = 0;

	template <class T>
	T* ExtractData()
	{
		T* p = dynamic_cast<T*>(this);
		if (p) return p; 
		else
		{
			if (m_pt) return m_pt->ExtractData<T>();
			else return 0;
		}
	}

protected:
	FEMaterialPoint*	m_pt;	//<! nested point data

public:
	static double time;	// time value
	static double dt; // time increment
};

//-----------------------------------------------------------------------------

class FEElasticMaterialPoint : public FEMaterialPoint
{
public:
	FEMaterialPoint* Copy()
	{
		FEElasticMaterialPoint* pt = new FEElasticMaterialPoint(*this);
		if (m_pt) pt->m_pt = m_pt->Copy();
		return pt;
	}

	void Serialize(Archive& ar)
	{
		if (ar.IsSaving())
		{
			ar << F << J << Q << s << avgJ << avgp;
		}
		else
		{
			ar >> F >> J >> Q >> s >> avgJ >> avgp;
		}

		if (m_pt) m_pt->Serialize(ar);
	}

	mat3ds Strain();

	mat3ds RightCauchyGreen();
	mat3ds LeftCauchyGreen();

	mat3ds pull_back(const mat3ds& A);
	mat3ds push_forward(const mat3ds& A);

public:
	void Init(bool bflag)
	{
		if (bflag)
		{
			F.unit();

			J = 1;

			s.zero();

			avgJ = avgp = 0;
		}

		if (m_pt) m_pt->Init(bflag);
	}

public:
	// position 
//	vec3d	r0;	//!< material position
//	vec3d	rt;	//!< spatial position

	// deformation data
	mat3d	F;	//!< deformation gradient
	double	J;			//!< determinant of F
	mat3d	Q;			//!< local material orientation

	// solid material data
	mat3ds		s;			//!< Cauchy stress

	// three-field element material data
	double		avgJ;		//!< average Jacobian (for three-field elements)
	double		avgp;		//!< average pressure (for three-field elements)
};

//-----------------------------------------------------------------------------

class FEPoroElasticMaterialPoint : public FEMaterialPoint
{
public:
	FEPoroElasticMaterialPoint(FEMaterialPoint* ppt) : FEMaterialPoint(ppt) {}

	FEMaterialPoint* Copy()
	{
		FEPoroElasticMaterialPoint* pt = new FEPoroElasticMaterialPoint(*this);
		if (m_pt) pt->m_pt = m_pt->Copy();
		return pt;
	}

	void Serialize(Archive& ar)
	{
		if (ar.IsSaving())
		{
			ar << m_p << m_gradp << m_w;
		}
		else
		{
			ar >> m_p >> m_gradp >> m_w;
		}

		if (m_pt) m_pt->Serialize(ar);
	}

	void Init(bool bflag)
	{
		if (bflag)
		{
			m_p = 0;
			m_gradp = vec3d(0,0,0);
			m_w = vec3d(0,0,0);
		}

		if (m_pt) m_pt->Init(bflag);
	}

public:
	// poro-elastic material data
	double		m_p;		//!< fluid pressure
	vec3d		m_gradp;	//!< spatial gradient of p
	vec3d		m_w;		//!< fluid flux
};

//-----------------------------------------------------------------------------

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

	void Serialize(Archive& ar)
	{
		if (ar.IsSaving())
		{
			ar << m_Sep;
			ar << (int) MAX_TERMS;
			for (int i=0; i<MAX_TERMS; ++i) ar << m_H[i] << m_Hp[i];
		}
		else
		{
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
