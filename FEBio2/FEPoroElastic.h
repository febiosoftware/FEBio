// FEPoroElastic.h: interface for the FEPoroElastic class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FEPOROELASTIC_H__C67341B3_B080_4E3B_8B34_D73AEF86BB33__INCLUDED_)
#define AFX_FEPOROELASTIC_H__C67341B3_B080_4E3B_8B34_D73AEF86BB33__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FECore/FEMaterial.h"
#include "FECore/FEElasticMaterial.h"

//-----------------------------------------------------------------------------
// Biphasic material point class.
//
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
			ar << m_p << m_gradp << m_w << m_pa << m_phiw;
		}
		else
		{
			ar >> m_p >> m_gradp >> m_w >> m_pa >> m_phiw;
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
			m_phiw = 1;
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
	double		m_phiw;		//!< porosity in current configuration
};

//-----------------------------------------------------------------------------
//! Base class for poroelastic materials.

class FEPoroElastic : public FENestedMaterial
{
public:
	FEPoroElastic();

	// returns a pointer to a new material point object
	FEMaterialPoint* CreateMaterialPointData() 
	{ 
		return new FEPoroElasticMaterialPoint(m_pBase->CreateMaterialPointData());
	}

public:
	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt);

	//! calculate tangent stiffness at material point
	virtual tens4ds Tangent(FEMaterialPoint& pt);

	//! calculate fluid flux
	virtual vec3d Flux(FEMaterialPoint& pt) = 0;

	//! permeability
	virtual void Permeability(double k[3][3], FEMaterialPoint& pt) = 0;

	//! tangent of permeability
	virtual tens4ds Tangent_Permeability(FEMaterialPoint& mp) = 0;

};

#endif // !defined(AFX_FEPOROELASTIC_H__C67341B3_B080_4E3B_8B34_D73AEF86BB33__INCLUDED_)
