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
#include "FEReactiveVEMaterialPoint.h"

//-----------------------------------------------------------------------------
//! This class implements a large deformation reactive viscoelastic material
//
class FEReactiveViscoelasticMaterial :	public FEElasticMaterial
{
public:
	//! default constructor
	FEReactiveViscoelasticMaterial(FEModel* pfem);
       
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
    
public:
    //! data initialization
    void Init();
    
	//! stress function
	mat3ds Stress(FEMaterialPoint& pt);
    
	//! tangent function
	tens4ds Tangent(FEMaterialPoint& pt);
    
    //! strain energy density function
    double StrainEnergyDensity(FEMaterialPoint& pt);
    
    //! cull generations
    void CullGenerations(FEMaterialPoint& pt);
    
    //! evaluate bond mass fraction for a given generation
    double BreakingBondMassFraction(FEMaterialPoint& pt, const int ig);
    
    //! evaluate bond mass fraction of reforming generation
    double ReformingBondMassFraction(FEMaterialPoint& pt);
    
    //! detect new generation
    bool NewGeneration(FEMaterialPoint& pt);
    
	//! returns a pointer to a new material point object
	FEMaterialPoint* CreateMaterialPointData();
    
private:
	FEPropertyT<FEElasticMaterial>	m_pBase;	//!< pointer to elastic solid material for polymer base
	FEPropertyT<FEElasticMaterial>	m_pBond;	//!< pointer to elastic solid material for reactive bonds
	FEPropertyT<FEBondRelaxation>   m_pRelx;    //!< pointer to bond relaxation material for reactive bonds
    
public:
    double	m_wmin;		//!< minimum value of relaxation
    int     m_btype;    //!< bond kinetics type
    int     m_ttype;    //!< bond breaking trigger type
    
    DECLARE_PARAMETER_LIST();
};

#endif /* defined(__FEBioMech__FEReactiveViscoelastic__) */
