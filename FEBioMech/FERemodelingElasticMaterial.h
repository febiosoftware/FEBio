#pragma once
#include "FEElasticMaterial.h"

//-----------------------------------------------------------------------------
//! Base class for solid supply.
//! These materials need to define the solid supply and tangent supply functions.
//! The solid supply has units of mass/(referential volume)/time
//!

class FESolidSupply : public FEMaterial
{
public:
	//! constructor
	FESolidSupply(FEModel* pfem) : FEMaterial(pfem) {}

	//! Initialization
	virtual void Init() {}
	
	//! solid supply
	virtual double Supply(FEMaterialPoint& pt) = 0;
	
	//! tangent of solute supply with respect to strain
	virtual mat3ds Tangent_Supply_Strain(FEMaterialPoint& mp) = 0;
	
	//! tangent of solute supply with respect to referential density
	virtual double Tangent_Supply_Density(FEMaterialPoint& mp) = 0;	
};

//-----------------------------------------------------------------------------
//! Material point data for remodeling elastic materials
class FERemodelingMaterialPoint : public FEMaterialPoint
{
public:
	FERemodelingMaterialPoint(FEMaterialPoint *pt) : FEMaterialPoint(pt) {}
    
	FEMaterialPoint* Copy();
    
	void Init(bool bflag);
    
	void Serialize(DumpFile& ar);

	void ShallowCopy(DumpStream& dmp, bool bsave);
    
public:
	double		m_sed;		//!< strain energy density
	double		m_dsed;		//!< derivative of strain energy density with mass density
	double		m_rhor;		//!< current referential mass density
	double		m_rhorp;	//!< referential mass density at previous time step
};

//-----------------------------------------------------------------------------
//! A material that wants to use the remodeling framework needs to implement 
//! this additional interface.
class FERemodelingInterface
{
public:
	//! calculate strain energy density at material point
	virtual double StrainEnergy(FEMaterialPoint& pt) = 0;

	//! calculate tangent of strain energy density with solid density at material point
	virtual double Tangent_SE_Density(FEMaterialPoint& pt) = 0;

	//! calculate tangent of stress with solid density at material point
	virtual mat3ds Tangent_Stress_Density(FEMaterialPoint& pt) = 0;
};

//-----------------------------------------------------------------------------
//! Material class for remodeling solids
class FERemodelingElasticMaterial : public FEElasticMaterial
{
public:
	//! default constructor
	FERemodelingElasticMaterial(FEModel* pfem) : FEElasticMaterial(pfem) {}
	
	//! strain energy density function
	double StrainEnergy(FEMaterialPoint& pt);
	
	//! stress function
	mat3ds Stress(FEMaterialPoint& pt);
	
	//! tangent function of stress with strain
	tens4ds Tangent(FEMaterialPoint& pt);
	
	//! tangent function of strain energy density with solid mass density
	double Tangent_SE_Density(FEMaterialPoint& pt);
	
	//! tangent function of stress with solid mass density
	mat3ds Tangent_Stress_Density(FEMaterialPoint& pt);
	
	//! data initialization and checking
	void Init();
	
	// returns a pointer to a new material point object
	FEMaterialPoint* CreateMaterialPointData()
	{
		return new FERemodelingMaterialPoint(m_pBase->CreateMaterialPointData());
	}
    
	// get the elastic material
	FEElasticMaterial* GetElasticMaterial() { return m_pBase; }

public:
	//! return number of material properties
	int Properties();

	//! return a material property
	FEMaterial* GetProperty(int i);

	//! find a material property index ( returns <0 for error)
	int FindPropertyIndex(const char* szname);

	//! set a material property (returns false on error)
	bool SetProperty(int i, FEMaterial* pm);
	
public:
	FEElasticMaterial*	m_pBase;		//!< pointer to elastic solid material
	FESolidSupply*		m_pSupp;		//!< pointer to solid supply material
	double				m_rhormin;		//!< minimum density
	double				m_rhormax;		//!< maximum density
	
public:
	// declare parameter list
	DECLARE_PARAMETER_LIST();
};
