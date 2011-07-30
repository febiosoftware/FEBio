#pragma once

#include "mat3d.h"
#include "DumpFile.h"
#include <vector>
using namespace std;

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

	virtual void Serialize(DumpFile& ar) = 0;

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
	FEElasticMaterialPoint()
	{
		F.zero();
		Q.unit();
		J = 1;
		s.zero();
		s0.zero();
	}

	FEMaterialPoint* Copy()
	{
		FEElasticMaterialPoint* pt = new FEElasticMaterialPoint(*this);
		if (m_pt) pt->m_pt = m_pt->Copy();
		return pt;
	}

	void Serialize(DumpFile& ar)
	{
		if (ar.IsSaving())
		{
			ar << F << J << Q << s << s0;
		}
		else
		{
			ar >> F >> J >> Q >> s >> s0;
		}

		if (m_pt) m_pt->Serialize(ar);
	}

	mat3ds Strain();

	mat3ds RightCauchyGreen();
	mat3ds LeftCauchyGreen ();

	mat3ds DevRightCauchyGreen();
	mat3ds DevLeftCauchyGreen ();

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
			s0.zero();

//			Q.unit();
		}

		if (m_pt) m_pt->Init(bflag);
	}

public:
	// position 
	vec3d	r0;	//!< material position
	vec3d	rt;	//!< spatial position

	// deformation data
	mat3d	F;	//!< deformation gradient
	double	J;			//!< determinant8 of F
	mat3d	Q;			//!< local material orientation

	// solid material data
	mat3ds		s;			//!< Cauchy stress
	mat3ds		s0;			//!< Initial stress (only used by linear solid solver)
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

	void Serialize(DumpFile& ar)
	{
		if (ar.IsSaving())
		{
			ar << m_p << m_gradp << m_w << m_pa;
		}
		else
		{
			ar >> m_p >> m_gradp >> m_w >> m_pa;
		}

		if (m_pt) m_pt->Serialize(ar);
	}

	void Init(bool bflag)
	{
		if (bflag)
		{
			m_p = m_pa = 0;
			m_gradp = vec3d(0,0,0);
			m_w = vec3d(0,0,0);
		}

		if (m_pt) m_pt->Init(bflag);
	}

public:
	// poro-elastic material data
	// The actual fluid pressure is the same as the effective fluid pressure
	// in a poroelastic material without solute(s).  The actual fluid pressure
	// is included here so that models that include both poroelastic and
	// solute-poroelastic domains produce plotfiles with consistent fluid
	// pressure fields.
	double		m_p;		//!< fluid pressure
	vec3d		m_gradp;	//!< spatial gradient of p
	vec3d		m_w;		//!< fluid flux
	double		m_pa;		//!< actual fluid pressure
};

//-----------------------------------------------------------------------------

class FEHeatMaterialPoint : public FEMaterialPoint
{
public:
	FEMaterialPoint* Copy()
	{
		FEHeatMaterialPoint* pt = new FEHeatMaterialPoint(*this);
		if (m_pt) pt->m_pt = m_pt->Copy();
		return pt;
	}

	void Serialize(DumpFile& ar)
	{
		if (m_pt) m_pt->Serialize(ar);
	}

	void Init(bool bflag)
	{
		if (m_pt) m_pt->Init(bflag);
	}

public:
	vec3d	m_q;	//!< heat flux
};

//-----------------------------------------------------------------------------

class FETrussMaterialPoint : public FEMaterialPoint
{
public:
	FEMaterialPoint* Copy()
	{
		FETrussMaterialPoint* pt = new FETrussMaterialPoint(*this);
		if (m_pt) pt->m_pt = m_pt->Copy();
		return pt;
	}

	void Serialize(DumpFile& ar)
	{
		if (m_pt) m_pt->Serialize(ar);
	}

	void Init(bool bflag)
	{
		if (m_pt) m_pt->Init(bflag);
		m_l = 1;
		m_tau = 0;
	}

public:
	double	m_l;	// strech
	double	m_tau;	// Kirchoff stress
};

//-----------------------------------------------------------------------------

class FEDiscreteMaterialPoint : public FEMaterialPoint
{
public:
	FEMaterialPoint* Copy()
	{
		FEDiscreteMaterialPoint* pt = new FEDiscreteMaterialPoint(*this);
		if (m_pt) pt->m_pt = m_pt->Copy();
		return pt;
	}

	void Serialize(DumpFile& ar)
	{
		if (m_pt) m_pt->Serialize(ar);
	}

	void Init(bool bflag)
	{
		if (m_pt) m_pt->Init(bflag);
	}
};

//-----------------------------------------------------------------------------

class FESolutePoroElasticMaterialPoint : public FEMaterialPoint
{
public:
	FESolutePoroElasticMaterialPoint(FEMaterialPoint* ppt) : FEMaterialPoint(ppt) {}
	
	FEMaterialPoint* Copy()
	{
		FESolutePoroElasticMaterialPoint* pt = new FESolutePoroElasticMaterialPoint(*this);
		if (m_pt) pt->m_pt = m_pt->Copy();
		return pt;
	}
	
	void Serialize(DumpFile& ar)
	{
		if (ar.IsSaving())
		{
			ar << m_p << m_gradp << m_w << m_pa << m_c << m_gradc << m_j << m_ca;
		}
		else
		{
			ar >> m_p >> m_gradp >> m_w >> m_pa >> m_c >> m_gradc >> m_j >> m_ca;
		}
		
		if (m_pt) m_pt->Serialize(ar);
	}
	
	void Init(bool bflag)
	{
		if (bflag)
		{
			m_p = m_pa = 0;
			m_gradp = vec3d(0,0,0);
			m_w = vec3d(0,0,0);
			m_c = m_ca = 0;
			m_gradc = vec3d(0,0,0);
			m_j = vec3d(0,0,0);
		}
		
		if (m_pt) m_pt->Init(bflag);
	}
	
public:
	// biphasic-solute material data
	double		m_p;		//!< effective fluid pressure
	vec3d		m_gradp;	//!< spatial gradient of p
	vec3d		m_w;		//!< fluid flux
	double		m_pa;		//!< actual fluid pressure
	// solute material data
	double		m_c;		//!< effective solute concentration
	vec3d		m_gradc;	//!< spatial gradient of c
	vec3d		m_j;		//!< solute molar flux
	double		m_ca;		//!< actual solute concentration
};

//-----------------------------------------------------------------------------
//! Material point data for mixtures
//!
class FEElasticMixtureMaterialPoint : public FEMaterialPoint
{
public:
	FEElasticMixtureMaterialPoint() { m_pt = new FEElasticMaterialPoint; }
	FEMaterialPoint* Copy()
	{
		FEElasticMixtureMaterialPoint* pt = new FEElasticMixtureMaterialPoint;
		pt->m_w = m_w;
		if (m_pt) pt->m_pt = m_pt->Copy();
		return pt;
	}

	void Init(bool bflag)
	{
		if (bflag)
		{
			for (int i=0; i<(int) m_w.size(); ++i) m_w[i] = 1.0;
		}
	}

	void Serialize(DumpFile& ar)
	{
		if (ar.IsSaving())
		{
			ar << m_w;
		}
		else
		{
			ar >> m_w;
		}
	}

public:
	vector<double>	m_w;	//!< material weights
};
