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
#include <string.h>
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
	FEMaterial();
	virtual ~FEMaterial() {}

	//! set/get material name
	void SetName(const char* sz) { strcpy(m_szname, sz); }
	const char* GetName() { return m_szname; }

	//! returns a pointer to a new material point object
	virtual FEMaterialPoint* CreateMaterialPointData() { return 0; };

	//! performs initialization and parameter checking
	virtual void Init(){}

	int GetID() { return m_nID; }

public:
	// this is the first GetParameterList function
	// so we have to declare it explicitly
	virtual FEParameterList* GetParameterList();

	// returns the type string of the material
	virtual const char* GetTypeString();

protected:
	char	m_szname[128];	//!< name of material
	int		m_nID;			//!< material ID
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
	FESolidMaterial*	m_pBase;	//!< pointer to base material

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
//! Base class for hydraulic permeability of porous materials.
//! These materials need to define the permeability and tangent permeability functions.
//!
class FEHydraulicPermeability : public FEMaterial
	{
	public:
		FEHydraulicPermeability() {m_phi0 = -1; }
		virtual ~FEHydraulicPermeability(){}
		
		//! hydraulic permeability
		virtual mat3ds Permeability(FEMaterialPoint& pt) = 0;
		
		//! tangent of hydraulic permeability with respect to strain
		virtual tens4ds Tangent_Permeability_Strain(FEMaterialPoint& mp) = 0;
		
		//! tangent of hydraulic permeability with respect to concentration
		mat3ds Tangent_Permeability_Concentration(FEMaterialPoint& mp);
		
		void Init();
		
	public:
		double	m_phi0;			//!< solid volume fraction in reference state
		
		// declare parameter list
		DECLARE_PARAMETER_LIST();
	};

//-----------------------------------------------------------------------------
//! Base class for biphasic materials.

class FEBiphasic : public FEMaterial
	{
	public:
		FEBiphasic();
		
		// returns a pointer to a new material point object
		FEMaterialPoint* CreateMaterialPointData() 
		{ 
			return new FEPoroElasticMaterialPoint(m_pSolid->CreateMaterialPointData());
		}
		
	public:
		void Init();
		
		//! calculate stress at material point
		mat3ds Stress(FEMaterialPoint& pt);
		
		//! calculate tangent stiffness at material point
		tens4ds Tangent(FEMaterialPoint& pt);
		
		//! calculate fluid flux
		vec3d Flux(FEMaterialPoint& pt);
		
		//! calculate actual fluid pressure
		double Pressure(FEMaterialPoint& pt);
		
		//! permeability
		void Permeability(double k[3][3], FEMaterialPoint& pt);
		
		//! tangent of permeability
		tens4ds Tangent_Permeability_Strain(FEMaterialPoint& pt);
		
		//! porosity
		double Porosity(FEMaterialPoint& pt);
		
		//! fluid density
		double FluidDensity() { return m_rhoTw; } 
		
	public:
		double						m_rhoTw;	//!< true fluid density
		FEElasticMaterial*			m_pSolid;	//!< pointer to elastic solid material
		FEHydraulicPermeability*	m_pPerm;	//!< pointer to permeability material
		
		// declare as registered
		DECLARE_REGISTERED(FEBiphasic);
		
		DECLARE_PARAMETER_LIST();
		
	};

#endif // !defined(AFX_FEMATERIAL_H__07F3E572_45B6_444E_A3ED_33FE9D18E82D__INCLUDED_)
