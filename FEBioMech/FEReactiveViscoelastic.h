//
//  FEReactiveViscoelastic.h
//  FEBioMech
//
//  Created by Gerard Ateshian on 8/25/14.
//  Copyright (c) 2014 febio.org. All rights reserved.
//

#ifndef __FEBioMech__FEReactiveViscoelastic__
#define __FEBioMech__FEReactiveViscoelastic__

#include "FEElasticMaterial.h"
#include "FEBondRelaxation.h"

//-----------------------------------------------------------------------------
//! Material point data for reactive viscoelastic materials
class FEReactiveVEMaterialPoint : public FEMaterialPoint
{
public:
	//! constructor
	FEReactiveVEMaterialPoint(FEMaterialPoint *pt) : FEMaterialPoint(pt) {}
    
	//! copy material point data
	FEMaterialPoint* Copy();
    
	//! Initialize material point data
	void Init(bool bflag);
    
	//! Serialize data to archive
	void Serialize(DumpFile& ar);
    
	//! data streaming
	void ShallowCopy(DumpStream& dmp, bool bsave);
    
public:
	// multigenerational material data
	vector <mat3d>  m_Fi;	//!< inverse of relative deformation gradient
	vector <double> m_Ji;	//!< determinant of Fi (store for efficiency)
	vector <double> m_tgen;	//!< generation time
    
};


//-----------------------------------------------------------------------------
//! This class implements a large deformation reactive viscoelastic material
//
class FEReactiveViscoelasticMaterial :	public FEElasticMaterial
{
public:
	//! default constructor
	FEReactiveViscoelasticMaterial(FEModel* pfem);
    
	//! Get a parameter
	FEParam* GetParameter(const ParamString& s);
    
	//! get the elastic base material
	FEElasticMaterial* GetBaseMaterial() { return m_pBase; }
    
	//! Set the base material
	void SetBaseMaterial(FEElasticMaterial* pbase) { m_pBase = pbase; }
    
	//! get the elastic bond material
	FEElasticMaterial* GetBondMaterial() { return m_pBond; }
    
	//! Set the base material
	void SetBondMaterial(FEElasticMaterial* pbond) { m_pBond = pbond; }
    
	//! Set the local coordinate system for a material point (overridden from FEMaterial)
	void SetLocalCoordinateSystem(FEElement& el, int n, FEMaterialPoint& mp);
    
	//! serialize data to/from dump file
	void Serialize(DumpFile& ar);
    
public:
	//! return number of properties
	int Properties();
    
	//! return a material property
	FECoreBase* GetProperty(int i);
    
	//! find a material property index ( returns <0 for error)
	virtual int FindPropertyIndex(const char* szname);
    
	//! set a material property (returns false on error)
	virtual bool SetProperty(int i, FECoreBase* pm);
    
public:
	//! data initialization
	void Init();
    
	//! stress function
	mat3ds Stress(FEMaterialPoint& pt);
    
	//! tangent function
	tens4ds Tangent(FEMaterialPoint& pt);
    
    //! strain energy density function
    double StrainEnergyDensity(FEMaterialPoint& pt);
    
	// returns a pointer to a new material point object
	FEMaterialPoint* CreateMaterialPointData();
    
private:
	FEElasticMaterial*	m_pBase;	//!< pointer to elastic solid material for polymer base
	FEElasticMaterial*	m_pBond;	//!< pointer to elastic solid material for reactive bonds
    FEBondRelaxation*   m_pRelx;    //!< pointer to bond relaxation material for reactive bonds
};

#endif /* defined(__FEBioMech__FEReactiveViscoelastic__) */
