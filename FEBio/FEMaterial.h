// FEMaterial.h: interface for the FEMaterial class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FEMATERIAL_H__07F3E572_45B6_444E_A3ED_33FE9D18E82D__INCLUDED_)
#define AFX_FEMATERIAL_H__07F3E572_45B6_444E_A3ED_33FE9D18E82D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FEElement.h"
#include "quatd.h"
#include "tens4d.h"
#include "LoadCurve.h"
#include "string.h"
#include "FE_enum.h"
#include "FECoordSysMap.h"
#include "FEMaterialFactory.h"
#include "FEParameterList.h"
#include "FEMaterialPoint.h"

#define INRANGE(x, a, b) ((x)>=(a) && (x)<=(b))

//-----------------------------------------------------------------------------
//! class to throw during the material initialization phase

class MaterialError
{
public:
	MaterialError(const char* sz, ...);

	const char* Error() { return m_szerr; }

protected:
	char	m_szerr[512];
};

//-----------------------------------------------------------------------------
//! Abstract base class for material types

//! From this class all other material classes are derived.

class FEMaterial  
{
public:
	FEMaterial() { m_szname[0] = 0; }
	virtual ~FEMaterial() {}

	//! set/get material name
	void SetName(const char* sz) { strcpy(m_szname, sz); }
	const char* GetName() { return m_szname; }

	//! returns a pointer to a new material point object
	virtual FEMaterialPoint* CreateMaterialPointData() { return 0; };

	//! performs initialization and parameter checking
	virtual void Init(){}

public:
	// this is the first GetParameterList function
	// so we have to declare it explicitly
	virtual FEParameterList* GetParameterList();

	// returns the type string of the material
	virtual const char* GetTypeString();

protected:
	char	m_szname[128];	//!< name of material
};

//-----------------------------------------------------------------------------
//! Base class for solid-materials.
//! These materials need to define the stress and tangent functions.
//!
class FESolidMaterial : public FEMaterial
{
public:
	//! calculate stress at material point
	virtual mat3ds Stress(FEMaterialPoint& pt) = 0;

	//! calculate tangent stiffness at material point
	virtual tens4ds Tangent(FEMaterialPoint& pt) = 0;

	//! return the bulk modulus
	virtual double BulkModulus() = 0;

	//! return the material density
	virtual double Density() = 0;
};

//-----------------------------------------------------------------------------
//! Base class for (hyper-)elastic materials

class FEElasticMaterial : public FESolidMaterial
{
public:
	FEElasticMaterial() { m_density = 1; m_pmap = 0; m_unstable = false;}
	~FEElasticMaterial(){ if(m_pmap) delete m_pmap; }

	virtual FEMaterialPoint* CreateMaterialPointData() { return new FEElasticMaterialPoint; }

	void Init();

	double Density() { return m_density; } 

public:
	double	m_density;	//!< material density
	bool	m_unstable;	//!< flag indicating whether material is unstable on its own

	FECoordSysMap*	m_pmap;	//!< local material coordinate system

	DECLARE_PARAMETER_LIST();
};

//-----------------------------------------------------------------------------
//! Base class for nested materials. A nested material describes whose material 
//! response depends on the formulation of another user specified material. For
//! instance, the FEPoroElastic and FEViscoElastic are two examples of nested
//! materials.

class FENestedMaterial : public FESolidMaterial
{
public:
	FENestedMaterial() { m_nBaseMat = -1; m_pBase = 0; }
	virtual ~FENestedMaterial(){}

	//! return solid component's bulk modulus
	double BulkModulus() { return m_pBase->BulkModulus(); }

	//! return solid component's density
	double Density () { return m_pBase->Density(); }

public:
	int					m_nBaseMat;	//!< material ID of base material (one-based!)
	FEElasticMaterial*	m_pBase;	//!< pointer to base material

	DECLARE_PARAMETER_LIST();
};

//-----------------------------------------------------------------------------
//! Base class for incompressible materials

//! In FEBio a decoupled strain energy function is defined for incompressible
//! materials and a slightly different stress calculation is therefor required.
//! All materials who use this formulation need to derive from this class
//!
//! \todo I'm not happy with this class. I need to figure out a different
//! way to identify a material as incompressible

class FEIncompressibleMaterial : public FEElasticMaterial
{
public:
	//! constructor
	/*! \param ntype the unique material ID */
	FEIncompressibleMaterial() { m_blaugon = false; m_atol = 0.01; }

	// dilatational component U and derivs
	// use these to obtain similar results than NIKE3D
//	double U  (double J) { return 0.5*m_K*log(J)*log(J); }
	double Up (double J) { return m_K*log(J)/J; }
	double Upp(double J) { return m_K*(1-log(J))/(J*J); }

	// use these for NIKE3D's Ogden material
//	double U  (double J) { return 0.25*m_K*(J*J - 2.0*log(J) - 1.0); }
//	double Up (double J) { return 0.5*m_K*(J - 1.0/J); }
//	double Upp(double J) { return 0.5*m_K*(1 + 1.0/(J*J)); }

	// Use these to obtain similar results than ABAQUS
//	double U  (double J) { return 0.5*m_K*(J-1)*(J-1); }
//	double Up (double J) { return m_K*(J-1); }
//	double Upp(double J) { return m_K; }

	// incompressibility constraint fnc and derivs
	double h  (double J) { return log(J); }
	double hp (double J) { return 1.0 / J; }
	double hpp(double J) { return -1.0 / (J*J); }

	//! return bulk modulus
	double BulkModulus() { return m_K; }

public:
	bool	m_blaugon;	//!< augmented lagrangian flag
	double	m_atol;		//!< augmented lagrangian tolerance
	double	m_K;		//!< bulk modulus

	DECLARE_PARAMETER_LIST();
};

//-----------------------------------------------------------------------------
// Base class for fiber materials.
//
class FEFiberMaterial : public FEMaterial
{
public:
	FEFiberMaterial()
	{
		m_plc = 0;
		m_lcna = -1;
		m_ascl = 0;

		m_c3 = m_c4 = m_c5 = 0;
		m_lam1 = 1;
	}

	//! Calculate the fiber stress
	mat3ds Stress(FEMaterialPoint& mp);

	//! Calculate the fiber tangent
	tens4ds Tangent(FEMaterialPoint& mp);

public:
	double	m_c3;		//!< Exponential stress coefficient
	double	m_c4;		//!< Fiber uncrimping coefficient
	double	m_c5;		//!< Modulus of straightened fibers

	double	m_lam1;		//!< fiber stretch for straightened fibers

	//--- time varying elastance active contraction data ---
	int		m_lcna;		//!< use active contraction or not
	double	m_ascl;		//!< activation scale factor
	double	m_ca0;		//!< intracellular calcium concentration
	double	m_beta;		//!< shape of peak isometric tension-sarcomere length relation
	double	m_l0;		//!< unloaded length
	double	m_refl;		//!< sarcomere length

	// we need load curve data for active contraction
	FELoadCurve* m_plc;	//!< pointer to current load curve values
};


//-----------------------------------------------------------------------------
//! Base class for transversely isotropic materials.

//! This class was created to simplify the implementation of the TransIso Mooney-
//! Rivlin and Veronda-Westmann models.
//! This class only stores some data

// TODO: don't derive it from FEIncompressibleMaterial. Materials that are both
// incompressible and have a fiber distribution should derive from 
// FEIncompressibleMaterial and from a fiber material class. Or perhaps they should
// have a fiber (or material axis) class as a member.

class FETransverselyIsotropic : public FEIncompressibleMaterial
{
public:
	//! constructor
	FETransverselyIsotropic() {}

public:
	FEFiberMaterial	m_fib;

	// declare parameter list
	DECLARE_PARAMETER_LIST();
};

//-----------------------------------------------------------------------------
#ifdef WIN32
inline double acosh(double x)
{
	if (x <= 1) return 0;
	return log(x + sqrt(x*x - 1));
}
#endif

#endif // !defined(AFX_FEMATERIAL_H__07F3E572_45B6_444E_A3ED_33FE9D18E82D__INCLUDED_)
