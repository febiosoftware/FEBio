// FERigid.h: interface for the FERigid class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FERIGID_H__42385DA7_ECE1_4862_B6E1_EFE5B4D4CC4B__INCLUDED_)
#define AFX_FERIGID_H__42385DA7_ECE1_4862_B6E1_EFE5B4D4CC4B__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FECore/FEMaterial.h"

//-----------------------------------------------------------------------------
//! Rigd body material data

//! Since rigid elements are skipped during the stiffness and residual calculations
//! we don't implement the Stress and Tangent functions

class FERigidMaterial : public FEElasticMaterial
{
public:
	FERigidMaterial()
	{
		m_com = 0;
		for (int i=0; i<6; ++i)
		{
			m_bc[i] =  0;	// rigid bodies are initially free
			m_fc[i] = -1;
			m_fs[i] =  0;
		}
		m_E = 1;
		m_v = 0;
		m_pmid = -1;
	}

public:
	double	m_E;		//!< Young's modulus
	double	m_v;		//!< Poisson's ratio

public:
	int		m_com;	//!< center of mass input flag
	vec3d	m_rc;	//!< center of mass

	int		m_bc[6];	//!< rigid boundary conditions

	int		m_fc[6];	//!< rigid force/torque loadcurve ID's
	double	m_fs[6];	//!< rigid force/torque scale factors

	int		m_nRB;		//!< rigid body ID
	int		m_pmid;		//!< parent material ID

public:
	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt){ return mat3ds(); }

	//! calculate tangent stiffness at material point
	virtual tens4ds Tangent(FEMaterialPoint& pt) { return tens4ds(); }

	//! data initialization
	void Init();

	//! serialization
	void Serialize(DumpFile& ar);

	//! return bulk modulus
	double BulkModulus() { return m_E/(3.0*(1.0 - 2.0*m_v));}

	// declare as registered
	DECLARE_REGISTERED(FERigidMaterial);

	// declare a parameter list
	DECLARE_PARAMETER_LIST();
};

#endif // !defined(AFX_FERIGID_H__42385DA7_ECE1_4862_B6E1_EFE5B4D4CC4B__INCLUDED_)
