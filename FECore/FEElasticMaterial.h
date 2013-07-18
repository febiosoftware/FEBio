#pragma once
#include "FEMaterial.h"
#include "FECoordSysMap.h"

//-----------------------------------------------------------------------------
// This class defines material point data for elastic materials.
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
	mat3ds SmallStrain();

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
	double	sed;		//!< strain energy density	\todo Is this a good place for this?
	double	rhor;		//!< current referential mass density
};

//-----------------------------------------------------------------------------
//! Base class for (hyper-)elastic materials

class FEElasticMaterial : public FESolidMaterial
{
public:
	FEElasticMaterial() { m_density = 1; m_molarmass = 0; m_pmap = 0; m_unstable = false;}
	~FEElasticMaterial(){ if(m_pmap) delete m_pmap; }

	virtual FEMaterialPoint* CreateMaterialPointData() { return new FEElasticMaterialPoint; }

	void Init();

	void Serialize(DumpFile& ar);

	//! Get the elastic component
	FEElasticMaterial* GetElasticMaterial() { return this; }

public:
	bool	m_unstable;	//!< flag indicating whether material is unstable on its own
	void	*pVoid;		//!< pointer to parent

	FECoordSysMap*	m_pmap;	//!< local material coordinate system

	DECLARE_PARAMETER_LIST();
};
